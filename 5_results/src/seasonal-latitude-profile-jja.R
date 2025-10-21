library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(Matrix)

source(Sys.getenv('UTILS_PARTIAL'))
source(Sys.getenv('DISPLAY_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--perturbations-augmented')
parser$add_argument('--samples-LNLGIS')
parser$add_argument('--samples-LNLGISSIF')
parser$add_argument('--xbase-monthly-2x25')
parser$add_argument('--output')
args <- parser$parse_args()

perturbations_base <- fst::read_fst(args$perturbations_augmented)
samples_LNLGIS <- readRDS(args$samples_LNLGIS)
samples_LNLGISSIF <- readRDS(args$samples_LNLGISSIF)
xbase_monthly_2x25 <- fst::read_fst(args$xbase_monthly_2x25)

stopifnot(
  ncol(samples_LNLGIS$alpha_df$value_samples) == ncol(samples_LNLGISSIF$alpha_df$value_samples)
)

# Seasonal profiles

xbase_seasonal_profiles <- xbase_monthly_2x25 %>%
  filter(between(latitude, -56, 80)) %>%
  mutate(
    month = lubridate::month(time),
    season = factor(case_when(
      between(month, 3, 5) ~ 'MAM',
      between(month, 6, 8) ~ 'JJA',
      between(month, 9, 11) ~ 'SON',
      TRUE ~ 'DJF'
    ), levels = c('DJF', 'MAM', 'JJA', 'SON'))
  ) %>%
  group_by(inventory, season, latitude) %>%
  summarise(
    # Average over the six-year study period
    value = GC_M2_DAY_TO_PGC_MONTH * sum(area * value) / 6,
    .groups = 'drop'
  ) %>%
  mutate(estimate = 'FLUXCOM')

perturbations_base_seasonal <- perturbations_base %>%
  mutate(
    month = lubridate::month(time),
    season = factor(case_when(
      between(month, 3, 5) ~ 'MAM',
      between(month, 6, 8) ~ 'JJA',
      between(month, 9, 11) ~ 'SON',
      TRUE ~ 'DJF'
    ), levels = c('DJF', 'MAM', 'JJA', 'SON')),
    inventory_season_latitude = interaction(
      inventory,
      season,
      latitude,
      drop = TRUE
    )
  )

perturbations_seasonal <- perturbations_base_seasonal %>%
  group_by(inventory_season_latitude, basis_vector) %>%
  summarise(
    # Average over the six-year study period
    value = KG_M2_S_TO_PGC_MONTH * sum(area * value) / 6,
    .groups = 'drop'
  ) %>%
  left_join(
    perturbations_base_seasonal %>%
      distinct(inventory_season_latitude, inventory, season, latitude),
    by = 'inventory_season_latitude'
  )

X_seasonal_profiles <- with(perturbations_seasonal, sparseMatrix(
  i = as.integer(inventory_season_latitude),
  j = as.integer(basis_vector),
  x = value,
  dims = c(nlevels(inventory_season_latitude), nlevels(basis_vector))
))

prior_emissions_seasonal <- perturbations_seasonal %>%
  group_by(inventory_season_latitude, inventory, season, latitude) %>%
  summarise(value = sum(value), .groups = 'drop') %>%
  select(-inventory_season_latitude) %>%
  mutate(estimate = 'Bottom-up')

posterior_emissions_seasonal_LNLGIS <- compute_posterior(
  prior_emissions_seasonal,
  X_seasonal_profiles,
  samples_LNLGIS,
  'v2.0 posterior'
)
posterior_emissions_seasonal_LNLGISSIF <- compute_posterior(
  prior_emissions_seasonal,
  X_seasonal_profiles,
  samples_LNLGISSIF,
  'v2.S posterior'
)

emissions_seasonal <- bind_rows(
  prior_emissions_seasonal,
  posterior_emissions_seasonal_LNLGIS,
  posterior_emissions_seasonal_LNLGISSIF
) %>%
  filter(inventory %in% c('bio_assim', 'bio_resp_tot')) %>%
  {
    x <- .

    bind_rows(
      x,
      x %>%
        group_by(estimate, season, latitude) %>%
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
  bind_rows(xbase_seasonal_profiles) %>%
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
    )
  )

output <- emissions_seasonal %>%
  filter(season == 'JJA') %>%
  ggplot(aes(x = latitude)) +
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
  facet_wrap(vars(inventory), scales = 'free_y', nrow = 1) +
  scale_x_continuous(breaks = seq(-30, 60, 30)) +
  scale_colour_manual(values = DISPLAY_SETTINGS$colour_key) +
  scale_fill_manual(values = DISPLAY_SETTINGS$colour_key) +
  scale_linetype_manual(values = DISPLAY_SETTINGS$linetype_key) +
  guides(fill = 'none') +
  labs(x = 'Latitude [°N]', y = 'Flux [PgC/month]', colour = NULL, fill = NULL, linetype = NULL) +
  theme(
    plot.margin = margin(t = 0, r = 0.1, b = 0, l = 0.1, unit = 'cm'),
    plot.title = element_blank(),
    axis.text.x = element_text(size = 8, colour = '#23373b'),
    axis.title.x = element_text(
      size = 10,
      colour = '#23373b',
      margin = margin(t = 0.2, r = 0, b = 0, l = 0, unit = 'cm')
    ),
    axis.text.y = element_text(size = 7, colour = '#23373b'),
    axis.title.y = element_text(size = 10, colour = '#23373b'),
    strip.text = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.position = 'bottom',
    legend.margin = margin(t = -0.2, r = 0, b = 0, l = 0, unit = 'cm')
  )

ggsave_base(
  args$output,
  output,
  width = DISPLAY_SETTINGS$supplement_full_width,
  height = 5.5
)
