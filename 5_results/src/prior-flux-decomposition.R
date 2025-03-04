library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(Matrix)
library(Rcpp)
library(patchwork)

rcpp_cache_dir <- Sys.getenv('RCPP_CACHE_DIR')
options(rcpp.cache.dir = if (rcpp_cache_dir == '') tempdir() else rcpp_cache_dir)

source(Sys.getenv('UTILS_PARTIAL'))
source(Sys.getenv('DISPLAY_PARTIAL'))
sourceCpp(Sys.getenv('UTILS_CPP_PARTIAL'))
sourceCpp(Sys.getenv('HMC_EXACT_CPP_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--prior')
parser$add_argument('--constraints')
parser$add_argument('--basis-vectors')
parser$add_argument('--perturbations-augmented')
parser$add_argument('--region')
parser$add_argument('--output')
args <- parser$parse_args()

prior <- readRDS(args$prior)
constraints <- readRDS(args$constraints)
basis_vectors <- fst::read_fst(args$basis_vectors)
perturbations_augmented <- fst::read_fst(args$perturbations_augmented)

perturbations <- perturbations_augmented %>%
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

n_all_alpha <- nrow(basis_vectors)
# NOTE(mgnb): values of alpha fixed to zero are given infinite precision
alpha_to_include <- is.finite(diag(prior$precision))

alpha_prior_mean <- prior$mean[alpha_to_include]
alpha_prior_precision <- prior$precision[alpha_to_include, alpha_to_include]

bio_regions <- sort(unique(basis_vectors$region))[1 : 11]
bio_inventories <- c('bio_assim', 'bio_resp_tot')

n_alpha <- length(alpha_prior_mean)
n_bio_regions <- length(bio_regions)

F_constraint <- rbind(
  constraints$F_sign,
  constraints$F_residual
)[, alpha_to_include]
g_constraint <- c(constraints$g_sign, constraints$g_residual)

get_emissions_parts <- function(filter_expr) {
  perturbations_region <- if (missing(filter_expr)) {
    perturbations
  } else {
    perturbations %>% filter({{ filter_expr }})
  }
  perturbations_region <- perturbations_region %>%
    group_by(inventory_minor_component_time, basis_vector) %>%
    summarise(
      value = KG_M2_S_TO_PGC_MONTH * sum(area * value),
      .groups = 'drop'
    ) %>%
    left_join(
      perturbations %>%
        distinct(inventory_minor_component_time, inventory, minor_component, time),
      by = 'inventory_minor_component_time'
    ) %>%
    mutate(inventory_minor_component_time = droplevels(inventory_minor_component_time))

  X_region <- with(perturbations_region, sparseMatrix(
    i = as.integer(inventory_minor_component_time),
    j = as.integer(basis_vector),
    x = value,
    dims = c(nlevels(inventory_minor_component_time), nlevels(basis_vector))
  ))[, alpha_to_include]

  prior_emissions <- perturbations_region %>%
    group_by(inventory_minor_component_time, inventory, minor_component, time) %>%
    summarise(value = sum(value), .groups = 'drop') %>%
    select(-inventory_minor_component_time) %>%
    mutate(estimate = 'Bottom-up')

  list(
    X_region = X_region,
    prior_emissions = prior_emissions
  )
}

sample_from_region <- function(
  region_name,
  mu,
  Q,
  n_samples
) {
  indices_subset <- with(
    basis_vectors[alpha_to_include, ],
    region == region_name
  )
  F_constraint_subset <- F_constraint[, indices_subset]
  keep_F <- rowSums(abs(F_constraint_subset)) != 0
  F_constraint_subset <- F_constraint_subset[keep_F, ]
  g_constraint_subset <- g_constraint[keep_F]

  mu_subset <- mu[indices_subset]
  Q_subset <- Q[indices_subset, indices_subset]

  list(
    indices = indices_subset,
    samples = sampleHmcConstrained(
      rep(0, ncol(Q_subset)),
      mu_subset,
      R = chol(Q_subset),
      F = F_constraint_subset,
      g = g_constraint_subset,
      totalTime = pi / 2,
      debug = FALSE,
      nSamples = n_samples
    )
  )
}

compute_emissions <- function(alpha_samples, filter_expr) {
  emissions_parts <- get_emissions_parts(filter_expr)
  X_region <- emissions_parts$X_region
  prior_emissions <- emissions_parts$prior_emissions

  prior_emissions_samples <- prior_emissions %>%
    mutate(
      estimate = 'Prior',
      value_prior = value,
      value = value_prior + as.vector(
        X_region %*% colMeans(alpha_samples)
      ),
      value_samples = value_prior + as.matrix(
        X_region %*% t(alpha_samples)
      ),
      value_q025 = matrixStats::rowQuantiles(value_samples, probs = 0.025),
      value_q975 = matrixStats::rowQuantiles(value_samples, probs = 0.975)
    ) %>%
    select(-value_prior)

  sample_ids <- sample.int(nrow(alpha_samples), 5)
  bind_rows(
    prior_emissions,
    prior_emissions_samples
  ) %>%
    {
      x <- .

      bind_rows(
        x %>%
          mutate(
            value_s001 = value_samples[, sample_ids[1]],
            value_s002 = value_samples[, sample_ids[2]],
            value_s003 = value_samples[, sample_ids[3]],
            value_s004 = value_samples[, sample_ids[4]],
            value_s005 = value_samples[, sample_ids[5]]
          ),
        x %>%
          filter(inventory %in% c('bio_assim', 'bio_resp_tot')) %>%
          group_by(estimate, time, minor_component) %>%
          summarise(
            value = sum(value),
            value_samples = t(colSums(value_samples)),
            .groups = 'drop'
          ) %>%
          mutate(
            inventory = 'nee',
            value_q025 = matrixStats::rowQuantiles(value_samples, probs = 0.025),
            value_q975 = matrixStats::rowQuantiles(value_samples, probs = 0.975),
            value_s001 = value_samples[, sample_ids[1]],
            value_s002 = value_samples[, sample_ids[2]],
            value_s003 = value_samples[, sample_ids[3]],
            value_s004 = value_samples[, sample_ids[4]],
            value_s005 = value_samples[, sample_ids[5]]
          )
      )
    } %>%
    filter(!(
      (inventory == 'ocean' & minor_component %in% c('linear', 'periodic') & estimate == 'Prior')
    )) %>%
    mutate(
      inventory = factor(c(
        'bio_assim' = 'GPP',
        'bio_resp_tot' = 'Respiration',
        'nee' = 'NEE',
        'ocean' = 'Ocean'
      )[inventory], levels = c(
        'GPP',
        'Respiration',
        'NEE',
        'Ocean'
      )),
      minor_component = factor(c(
        'linear' = 'Linear',
        'periodic' = 'Seasonal',
        'residual' = 'Residual'
      )[as.character(minor_component)], levels = c(
        'Linear',
        'Seasonal',
        'Residual'
      ))
    )
}

plot_emissions <- function(alpha_samples, filter_expr) {
  emissions <- compute_emissions(alpha_samples, filter_expr)
  wrap_plots(lapply(sort(unique(emissions$inventory)), function(inventory_i) {
    ggplot(
      emissions %>%
        filter(inventory == inventory_i),
      aes(time)
    ) +
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
          y = value_s001,
          colour = estimate
        ),
        linetype = 'solid',
        linewidth = 0.4,
        alpha = 0.4
      ) +
      geom_line(
        mapping = aes(
          y = value_s002,
          colour = estimate
        ),
        linetype = 'solid',
        linewidth = 0.4,
        alpha = 0.4
      ) +
      geom_line(
        mapping = aes(
          y = value_s003,
          colour = estimate
        ),
        linetype = 'solid',
        linewidth = 0.4,
        alpha = 0.4
      ) +
      geom_line(
        mapping = aes(
          y = value_s004,
          colour = estimate
        ),
        linetype = 'solid',
        linewidth = 0.4,
        alpha = 0.4
      ) +
      geom_line(
        mapping = aes(
          y = value_s005,
          colour = estimate
        ),
        linetype = 'solid',
        linewidth = 0.4,
        alpha = 0.4
      ) +
      geom_line(
        mapping = aes(
          y = value,
          colour = estimate,
          linetype = estimate
        ),
        linewidth = 0.4
      ) +
      facet_grid(minor_component ~ ., scales = 'free_y') +
      scale_x_date(date_labels = '%Y-%m') +
      scale_colour_manual(values = DISPLAY_SETTINGS$colour_key) +
      scale_fill_manual(values = DISPLAY_SETTINGS$colour_key) +
      scale_linetype_manual(values = DISPLAY_SETTINGS$linetype_key) +
      labs(x = 'Time', y = 'Flux [PgC/month]', colour = NULL, fill = NULL, linetype = NULL) +
      guides(fill = 'none', linetype = 'none') +
      ggtitle(inventory_i)
  }), ncol = 2, nrow = 2, guides = 'collect') &
    theme(
      plot.margin = margin(t = 0.1, r = 0.1, b = 0, l = 0.1, unit = 'cm'),
      plot.title = element_text(
        size = 10,
        margin = margin(0, 0, 5.5, 0, unit = 'points')
      ),
      strip.text = element_text(size = 8),
      axis.text.x = element_text(size = 8, colour = '#23373b'),
      axis.title.x = element_text(
        size = 10,
        colour = '#23373b',
        margin = margin(t = 0.2, r = 0, b = 0, l = 0, unit = 'cm')
      ),
      axis.text.y = element_text(size = 7, colour = '#23373b'),
      axis.title.y = element_text(size = 10, colour = '#23373b'),
      legend.position = 'bottom',
      legend.margin = margin(t = -0.2, r = 0, b = 0, l = 0, unit = 'cm')
    )
}

n_prior_samples <- 1024
alpha_samples <- matrix(0, n_prior_samples, sum(alpha_to_include))
for (region_name in levels(basis_vectors$region)) {
  region_samples <- sample_from_region(
    region_name,
    alpha_prior_mean,
    alpha_prior_precision,
    n_prior_samples
  )
  alpha_samples[, region_samples$indices] <- region_samples$samples
}

output <- if (args$region == 'global') {
  plot_emissions(alpha_samples)
} else {
  plot_emissions(alpha_samples, region == args$region)
}

ggsave_base(
  args$output,
  output,
  width = DISPLAY_SETTINGS$supplement_full_width,
  height = 14.5
)
