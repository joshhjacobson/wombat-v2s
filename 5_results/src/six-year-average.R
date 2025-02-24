library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(Matrix)

source(Sys.getenv('UTILS_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--perturbations-augmented')
parser$add_argument('--samples-LNLGIS')
parser$add_argument('--samples-LNLGISSIF')
parser$add_argument('--xbase-monthly-2x25')
parser$add_argument('--output')
args <- parser$parse_args()

perturbations_augmented <- fst::read_fst(args$perturbations_augmented)
samples_LNLGIS <- readRDS(args$samples_LNLGIS)
samples_LNLGISSIF <- readRDS(args$samples_LNLGISSIF)

six_year_average_xbase <- fst::read_fst(args$xbase_monthly_2x25) %>%
  group_by(inventory, longitude, latitude, time) %>%
  summarise(
    value = GC_DAY_TO_KGCO2_YEAR * sum(value),
    .groups = 'drop'
  ) %>%
  group_by(inventory, longitude, latitude) %>%
  summarise(
    value = mean(value), .groups = 'drop'
  ) %>%
  mutate(
    estimate = 'FLUXCOM'
  )

perturbations_augmented <- perturbations_augmented %>%
  filter(
    inventory %in% c('bio_assim', 'bio_resp_tot')
  ) %>%
  mutate(
    value = PER_SECONDS_TO_PER_YEAR * value
  )

inventory_location <- perturbations_augmented %>%
  distinct(inventory, longitude, latitude) %>%
  arrange(inventory, longitude, latitude) %>%
  mutate(inventory_location_index = 1 : n())

X_alpha_to_six_year_mean <- with(
  perturbations_augmented %>%
    group_by(inventory, longitude, latitude, basis_vector) %>%
    summarise(value = sum(value) / 72, .groups = 'drop') %>%
    left_join(inventory_location, by = c('inventory', 'longitude', 'latitude')),
  sparseMatrix(
    i = inventory_location_index,
    j = as.integer(basis_vector),
    x = value,
    dims = c(nrow(inventory_location), nlevels(basis_vector))
  )
)

six_year_average_bottom_up <- inventory_location %>%
  select(-inventory_location_index) %>%
  left_join(
    perturbations_augmented %>%
      group_by(inventory, longitude, latitude, time) %>%
      summarise(
        value = sum(value), .groups = 'drop'
      ) %>%
      group_by(inventory, longitude, latitude) %>%
      summarise(
        value = mean(value), .groups = 'drop'
      ),
    by = c('inventory', 'longitude', 'latitude')
  )

compute_tilde <- function(samples) {
  design_matrix <- X_alpha_to_six_year_mean[, as.integer(samples$alpha_df$basis_vector)]
  alpha_mean <- samples$alpha_df$value
  alpha_samples <- samples$alpha_df$value_samples[
    ,
    floor(seq(1, ncol(samples$alpha_df$value_samples), length.out = 500))
  ]
  list(
    value_tilde_mean = as.vector(design_matrix %*% alpha_mean),
    value_tilde_samples = as.matrix(design_matrix %*% alpha_samples)
  )
}

tilde_LNLGIS <- compute_tilde(samples_LNLGIS)
tilde_LNLGISSIF <- compute_tilde(samples_LNLGISSIF)
tilde_difference <- list(
  value_tilde_mean = tilde_LNLGISSIF$value_tilde_mean - tilde_LNLGIS$value_tilde_mean,
  value_tilde_samples = tilde_LNLGISSIF$value_tilde_samples - tilde_LNLGIS$value_tilde_samples
)

output <- bind_rows(
  six_year_average_bottom_up %>%
    mutate(estimate = 'Bottom-up'),
  six_year_average_bottom_up %>%
    mutate(
      estimate = 'LNLGIS',
      value_samples = value + tilde_LNLGIS$value_tilde_samples,
      value = value + tilde_LNLGIS$value_tilde_mean,
      value_scale = matrixStats::rowSds(value_samples)
    ),
  six_year_average_bottom_up %>%
    mutate(
      estimate = 'LNLGISSIF',
      value_samples = value + tilde_LNLGISSIF$value_tilde_samples,
      value = value + tilde_LNLGISSIF$value_tilde_mean,
      value_scale = matrixStats::rowSds(value_samples)
    ),
  six_year_average_bottom_up %>%
    mutate(
      estimate = 'WOMBAT Difference',
      value_samples = tilde_difference$value_tilde_samples,
      value = tilde_difference$value_tilde_mean,
      value_scale = matrixStats::rowSds(value_samples)
    )
) %>%
  {
    x <- .

    bind_rows(
      x,
      x %>%
        group_by(estimate, longitude, latitude) %>%
        summarise(
          value_samples = t(colSums(value_samples)),
          value = sum(value),
          .groups = 'drop'
        ) %>%
        mutate(
          inventory = 'nee',
          value_scale = matrixStats::rowSds(value_samples)
        )
    )
  } %>%
  {
    x <- .

    bind_rows(
      x,
      x %>%
        filter(estimate == 'LNLGISSIF') %>%
        left_join(
          six_year_average_xbase %>%
            select(inventory, longitude, latitude, value_xbase = value),
          by = c('inventory', 'longitude', 'latitude')
        ) %>%
        mutate(
          estimate = 'FLUXCOM Difference',
          value = value - value_xbase
        ) %>%
        select(-c(value_xbase, value_scale))
    )
  } %>%
  select(-value_samples) %>%
  bind_rows(six_year_average_xbase)

fst::write_fst(output, args$output)
