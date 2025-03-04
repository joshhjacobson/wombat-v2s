library(argparse)
library(Matrix)
library(dplyr, warn.conflicts = FALSE)

source(Sys.getenv('UTILS_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--true-alpha')
parser$add_argument('--samples')
parser$add_argument('--perturbations-augmented')
parser$add_argument('--output')
args <- parser$parse_known_args()[[1]]

samples <- readRDS(args$samples)
perturbations_augmented <- fst::read_fst(args$perturbations_augmented)
true_alpha <- if (!is.null(args$true_alpha)) fst::read_fst(args$true_alpha) else NULL

perturbations_augmented <- perturbations_augmented %>%
  mutate(
    inventory_region_time = interaction(
      inventory,
      region,
      time,
      drop = TRUE
    )
  )

flux_aggregators <- perturbations_augmented %>%
  group_by(inventory_region_time, basis_vector) %>%
  summarise(
    value = KG_M2_S_TO_PGC_YEAR * sum(area * value),
    .groups = 'drop'
  ) %>%
  left_join(
    perturbations_augmented %>%
      distinct(inventory_region_time, inventory, region, time),
    by = 'inventory_region_time'
  )

X_aggregators <- with(flux_aggregators, sparseMatrix(
  i = as.integer(inventory_region_time),
  j = as.integer(basis_vector),
  x = value,
  dims = c(nlevels(inventory_region_time), nlevels(basis_vector))
))

prior_emissions <- flux_aggregators %>%
  group_by(inventory_region_time, inventory, region, time) %>%
  summarise(flux_mean = sum(value), .groups = 'drop') %>%
  select(-inventory_region_time) %>%
  mutate(estimate = 'Bottom-up')

true_emissions <- prior_emissions %>%
  mutate(estimate = 'Truth')

if (!is.null(true_alpha)) {
  true_emissions <- true_emissions %>%
    mutate(
      flux_mean = flux_mean + as.vector(
        X_aggregators[, as.integer(true_alpha$basis_vector)]
        %*% true_alpha$value
      )
    )
}

posterior_emissions <- prior_emissions %>%
  mutate(
    estimate = 'Posterior',
    flux_prior = flux_mean,
    flux_mean = flux_prior + as.vector(
      X_aggregators[, as.integer(samples$alpha_df$basis_vector)]
      %*% samples$alpha_df$value
    ),
    flux_samples = flux_prior + as.matrix(
      X_aggregators[, as.integer(samples$alpha_df$basis_vector)]
      %*% samples$alpha_df$value_samples
    )
  ) %>%
  select(-flux_prior)

flux_aggregates_samples <- bind_rows(
  true_emissions,
  prior_emissions,
  posterior_emissions
) %>%
  {
    x <- .

    bind_rows(
      x,
      x %>%
        filter(inventory %in% c('bio_assim', 'bio_resp_tot')) %>%
        group_by(estimate, region, time) %>%
        summarise(
          flux_mean = sum(flux_mean),
          flux_samples = t(colSums(flux_samples)),
          .groups = 'drop'
        ) %>%
        mutate(inventory = 'nee'),
      x %>%
        group_by(estimate, region, time) %>%
        summarise(
          flux_mean = sum(flux_mean),
          flux_samples = t(colSums(flux_samples)),
          .groups = 'drop'
        ) %>%
        mutate(inventory = 'total')
    )
  } %>%
  mutate(
    inventory = factor(c(
      'bio_assim' = 'GPP',
      'bio_resp_tot' = 'Respiration',
      'nee' = 'NEE',
      'ocean' = 'Ocean',
      'total' = 'Total'
    )[inventory], levels = c(
      'GPP',
      'Respiration',
      'NEE',
      'Ocean',
      'Total'
    ))
  )

saveRDS(flux_aggregates_samples, args$output)
