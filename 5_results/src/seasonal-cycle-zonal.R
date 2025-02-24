library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(Matrix)

source(Sys.getenv('UTILS_PARTIAL'))
source(Sys.getenv('DISPLAY_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--perturbations-augmented-zonal')
parser$add_argument('--samples-LNLGIS')
parser$add_argument('--samples-LNLGISSIF')
parser$add_argument('--xbase-monthly-2x25-zonal')
parser$add_argument('--output')
args <- parser$parse_args()

perturbations_zonal <- fst::read_fst(args$perturbations_augmented_zonal)
samples_LNLGIS <- readRDS(args$samples_LNLGIS)
samples_LNLGISSIF <- readRDS(args$samples_LNLGISSIF)
xbase_monthly_2x25_zonal <- fst::read_fst(args$xbase_monthly_2x25_zonal)

stopifnot(
  ncol(samples_LNLGIS$alpha_df$value_samples) == ncol(samples_LNLGISSIF$alpha_df$value_samples)
)

perturbations_zonal <- perturbations_zonal %>%
  mutate(
    month = lubridate::month(time),
    inventory_month = interaction(inventory, month, drop = TRUE)
  )

xbase_monthly_2x25_zonal <- xbase_monthly_2x25_zonal %>%
  mutate(
    month = lubridate::month(time)
  )

zones <- within(REGION_PLOT_SETTINGS, rm('global'))
emissions_zonal <- lapply(zones, function(zonal_band) {
  xbase_zone <- xbase_monthly_2x25_zonal %>%
    filter(
      latitude_bottom >= zonal_band$latitude_lower,
      latitude_bottom < zonal_band$latitude_upper
    ) %>%
    group_by(inventory, month) %>%
    summarise(
      # Average over the six-year study period
      value = GC_M2_DAY_TO_PGC_MONTH * sum(area * value) / 6,
      .groups = 'drop'
    ) %>%
    mutate(estimate = 'FLUXCOM')

  perturbations_zone <- perturbations_zonal %>%
    filter(
      latitude_bottom >= zonal_band$latitude_lower,
      latitude_bottom < zonal_band$latitude_upper
    ) %>%
    group_by(inventory_month, basis_vector) %>%
    summarise(
      # Average over the six-year study period
      value = KG_M2_S_TO_PGC_MONTH * sum(area * value) / 6,
      .groups = 'drop'
    ) %>%
    left_join(
      perturbations_zonal %>%
        distinct(inventory_month, inventory, month),
      by = 'inventory_month'
    )

  X_zone <- with(perturbations_zone, sparseMatrix(
    i = as.integer(inventory_month),
    j = as.integer(basis_vector),
    x = value,
    dims = c(nlevels(inventory_month), nlevels(basis_vector))
  ))

  prior_emissions <- perturbations_zone %>%
    group_by(inventory_month, inventory, month) %>%
    summarise(value = sum(value), .groups = 'drop') %>%
    select(-inventory_month) %>%
    mutate(estimate = 'Bottom-up')

  posterior_emissions_LNLGIS <- compute_posterior(prior_emissions, X_zone, samples_LNLGIS, 'v2.0 posterior')
  posterior_emissions_LNLGISSIF <- compute_posterior(prior_emissions, X_zone, samples_LNLGISSIF, 'v2.S posterior')

  emissions_zone <- bind_rows(
    prior_emissions,
    posterior_emissions_LNLGIS,
    posterior_emissions_LNLGISSIF
  ) %>%
    filter(inventory %in% c('bio_assim', 'bio_resp_tot')) %>%
    {
      x <- .

      bind_rows(
        x,
        x %>%
          group_by(estimate, month) %>%
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
    bind_rows(xbase_zone) %>%
    mutate(zone = zonal_band$numeric_title)
}) %>%
  bind_rows() %>%
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
    estimate = factor(
      estimate,
      levels = c('FLUXCOM', 'Bottom-up', 'v2.0 posterior', 'v2.S posterior')
    ),
    zone = factor(
      zone,
      levels = sapply(zones, function(zonal_band) zonal_band$numeric_title)
    )
  )

output <- emissions_zonal %>%
  ggplot(aes(x = month)) +
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
    linewidth = 0.5
  ) +
  ggh4x::facet_grid2(zone ~ inventory, scales = 'free_y', independent = 'y') +
  scale_y_continuous(n.breaks = 4) +
  scale_x_continuous(
    breaks = seq(1, 11, 2),
    labels = month.abb[c(TRUE, FALSE)]
  ) +
  scale_colour_manual(values = DISPLAY_SETTINGS$colour_key) +
  scale_fill_manual(values = DISPLAY_SETTINGS$colour_key) +
  scale_linetype_manual(values = DISPLAY_SETTINGS$linetype_key) +
  guides(fill = 'none') +
  labs(x = 'Month', y = 'Flux [PgC/month]', colour = NULL, fill = NULL, linetype = NULL) +
  theme(
    plot.margin = margin(t = 0, r = 0.1, b = 0, l = 0.1, unit = 'cm'),
    plot.title = element_blank(),
    axis.text.x = element_text(size = 9, colour = '#23373b'),
    axis.title.x = element_text(
      size = 10,
      colour = '#23373b',
      margin = margin(t = 0.2, r = 0, b = 0, l = 0, unit = 'cm')
    ),
    axis.text.y = element_text(size = 7, colour = '#23373b'),
    axis.title.y = element_text(size = 10, colour = '#23373b'),
    strip.text.x = element_text(size = 11),
    strip.text.y = element_text(size = 9),
    legend.text = element_text(size = 10),
    legend.position = 'bottom',
    legend.margin = margin(t = -0.2, r = 0, b = 0, l = 0, unit = 'cm')
  )

ggsave_base(
  args$output,
  output,
  width = DISPLAY_SETTINGS$supplement_full_width,
  height = 18
)
