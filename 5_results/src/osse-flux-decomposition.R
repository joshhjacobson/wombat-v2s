library(argparse)
library(Matrix)
library(dplyr, warn.conflicts = FALSE)

source(Sys.getenv('UTILS_PARTIAL'))
source(Sys.getenv('DISPLAY_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--perturbations-augmented')
parser$add_argument('--samples-fixresp-wsif')
parser$add_argument('--samples-fixresp-wosif')
parser$add_argument('--samples-freeresp-wsif')
parser$add_argument('--samples-freeresp-wosif')
parser$add_argument('--true-alpha')
parser$add_argument('--region')
parser$add_argument('--trim-labels', type = 'logical', default = FALSE)
parser$add_argument('--paper', type = 'logical', default = FALSE)
parser$add_argument('--output')
args <- parser$parse_known_args()[[1]]


true_alpha <- if (!is.null(args$true_alpha)) fst::read_fst(args$true_alpha) else NULL
perturbations <- fst::read_fst(args$perturbations_augmented)

if (args$region != 'global') {
  perturbations <- perturbations %>%
    filter(region == args$region)
}

perturbations <- perturbations %>%
  filter(inventory %in% c('bio_assim', 'bio_resp_tot')) %>%
  mutate(
    minor_component = factor(
      case_when(
        component == 'residual' ~ 'residual',
        component %in% c('intercept', 'trend') ~ 'linear',
        TRUE ~ 'periodic'
      ),
      c('linear', 'periodic', 'residual')
    ),
    inventory_minor_component_time = interaction(
      inventory,
      minor_component,
      time,
      drop = TRUE
    )
  )

perturbations_region <- perturbations %>%
  group_by(inventory_minor_component_time, basis_vector) %>%
  summarise(
    value = KG_M2_S_TO_PGC_MONTH * sum(area * value),
    .groups = 'drop'
  ) %>%
  left_join(
    perturbations %>%
      distinct(inventory_minor_component_time, inventory, minor_component, time),
    by = 'inventory_minor_component_time'
  )

X_region <- with(perturbations_region, sparseMatrix(
  i = as.integer(inventory_minor_component_time),
  j = as.integer(basis_vector),
  x = value,
  dims = c(nlevels(inventory_minor_component_time), nlevels(basis_vector))
))

prior_emissions <- perturbations_region %>%
  group_by(inventory_minor_component_time, inventory, minor_component, time) %>%
  summarise(value = sum(value), .groups = 'drop') %>%
  select(-inventory_minor_component_time) %>%
  mutate(estimate = 'Bottom-up')

true_emissions <- prior_emissions %>%
  mutate(estimate = 'Truth')

if (!is.null(true_alpha)) {
  log_debug('Adjusting bottom-up with true alpha from {args$true_alpha}')
  true_emissions <- true_emissions %>%
    mutate(
      value = value + as.vector(
        X_region[, as.integer(true_alpha$basis_vector)]
        %*% true_alpha$value
      )
    )
}

list_samples <- list(
  list(name = 'Without SIF, fixed RLT', path = args$samples_fixresp_wosif),
  list(name = 'With SIF, fixed RLT', path = args$samples_fixresp_wsif),
  list(name = 'Without SIF, inferred RLT', path = args$samples_freeresp_wosif),
  list(name = 'With SIF, inferred RLT', path = args$samples_freeresp_wsif)
)
list_samples <- list_samples[sapply(list_samples, function(x) !is.null(x$path))]
if (length(list_samples) < 1 || length(list_samples) > 2) {
  stop('Must provide at least one set of samples, and not more than two for comparison.')
}

posterior_emissions <- lapply(list_samples, function(samples_i) {
  if (is.null(samples_i$path)) return(NULL)
  samples <- readRDS(samples_i$path)
  compute_posterior(prior_emissions, X_region, samples, samples_i$name)
}) %>% bind_rows()

sample_names <- sapply(list_samples, function(x) x$name)
if (args$trim_labels) {
  sample_names <- sapply(sample_names, function(x) {
    strsplit(x, ',', fixed = TRUE)[[1]][1]
  }) %>% unique()

  posterior_emissions <- posterior_emissions %>%
    mutate(
      estimate = sapply(strsplit(estimate, ',', fixed = TRUE), function(x) x[1])
    )
}

emissions <- bind_rows(
  prior_emissions,
  posterior_emissions,
  true_emissions
) %>%
  {
    x <- .

    bind_rows(
      x,
      x %>%
        group_by(estimate, time, minor_component) %>%
        summarise(
          value = sum(value),
          value_samples = t(colSums(value_samples)),
          .groups = 'drop'
        ) %>%
        mutate(
          inventory = 'nee',
          value_q025 = matrixStats::rowQuantiles(value_samples, probs = 0.025),
          value_q975 = matrixStats::rowQuantiles(value_samples, probs = 0.975)
        )
    )
  } %>%
  filter(!(
    inventory == 'bio_resp_tot' & minor_component == 'linear' & grepl('fixed RLT', estimate, fixed = TRUE)
  )) %>%
  mutate(
    inventory = factor(c(
      'bio_assim' = 'GPP',
      'bio_resp_tot' = 'Respiration',
      'nee' = 'NEE'
    )[inventory], levels = c(
      'GPP',
      'Respiration',
      'NEE'
    )),
    minor_component = factor(c(
      'linear' = 'Linear',
      'periodic' = 'Seasonal',
      'residual' = 'Residual'
    )[as.character(minor_component)], levels = c(
      'Linear',
      'Seasonal',
      'Residual'
    )),
    estimate = factor(
      estimate,
      levels = c(
        'Truth',
        'Bottom-up',
        sample_names
      )
    )
  )

output <- ggplot(emissions, aes(x = time)) +
  geom_ribbon(
    mapping = aes(
      ymin = value_q025,
      ymax = value_q975,
      fill = estimate
    ),
    alpha = 0.3
  ) +
  geom_line(
    mapping = aes(
      y = value,
      colour = estimate,
      linetype = estimate
    ),
    linewidth = 0.5,
  ) +
  ggh4x::facet_grid2(minor_component ~ inventory, scales = 'free_y', independent = 'y') +
  scale_x_date(date_labels = '%Y-%m') +
  scale_colour_manual(values = DISPLAY_SETTINGS$colour_key) +
  scale_fill_manual(values = DISPLAY_SETTINGS$colour_key) +
  scale_linetype_manual(values = DISPLAY_SETTINGS$linetype_key) +
  labs(
    x = 'Time',
    y = 'Flux [PgC/month]',
    title = sprintf('Decomposition of %s fluxes', args$region),
    colour = NULL,
    fill = NULL,
    linetype = NULL
  ) +
  guides(fill = 'none') +
  theme(
    plot.margin = margin(t = 0, r = 0.1, b = 0, l = 0.1, unit = 'cm'),
    plot.title = if (args$paper) {
      element_blank()
    } else {
      element_text(
        hjust = 0.5,
        size = 11,
        margin = margin(0.1, 0, 0, 0, unit = 'cm')
      )
    },
    axis.text.x = element_text(size = 8, colour = '#23373b'),
    axis.title.x = element_text(
      size = 10,
      colour = '#23373b',
      margin = margin(t = 0.2, r = 0, b = 0, l = 0, unit = 'cm')
    ),
    axis.text.y = element_text(size = 7, colour = '#23373b'),
    axis.title.y = element_text(size = 10, colour = '#23373b'),
    strip.text.x = element_text(size = 10),
    strip.text.y = element_text(size = 9),
    legend.text = element_text(size = 10),
    legend.position = 'bottom',
    legend.margin = margin(t = -0.2, r = 0, b = 0, l = 0, unit = 'cm')
  )

ggsave_base(
  args$output,
  output,
  width = DISPLAY_SETTINGS$supplement_full_width,
  height = 8.5
)
