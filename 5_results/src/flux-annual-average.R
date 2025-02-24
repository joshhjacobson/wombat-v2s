library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(lubridate, warn.conflicts = FALSE)
library(Matrix)

source(Sys.getenv('UTILS_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--perturbations-augmented')
parser$add_argument('--samples-LNLGIS')
parser$add_argument('--samples-LNLGISSIF')
parser$add_argument('--xbase-monthly-2x25')
parser$add_argument('--output')
args <- parser$parse_args()

compute_posterior_mean_sd <- function(prior, design_matrix, samples, posterior_name) {
  prior %>%
    mutate(
      estimate = posterior_name,
      value_prior = value,
      value = value_prior + as.vector(
        design_matrix[, as.integer(samples$alpha_df$basis_vector)]
        %*% samples$alpha_df$value
      ),
      value_samples = value_prior + as.matrix(
        design_matrix[, as.integer(samples$alpha_df$basis_vector)]
        %*% samples$alpha_df$value_samples
      ),
      value_sd = matrixStats::rowSds(value_samples)
    ) %>%
    select(-value_prior)
}


perturbations_base <- fst::read_fst(args$perturbations_augmented)
samples_LNLGIS <- readRDS(args$samples_LNLGIS)
samples_LNLGISSIF <- readRDS(args$samples_LNLGISSIF)
xbase_monthly_2x25 <- fst::read_fst(args$xbase_monthly_2x25)

stopifnot(
  ncol(samples_LNLGIS$alpha_df$value_samples) == ncol(samples_LNLGISSIF$alpha_df$value_samples)
)

xbase_annual_average <- xbase_monthly_2x25 %>%
  mutate(year = year(time)) %>%
  group_by(inventory, year) %>%
  summarise(
    # sum over months gives PgC/year
    value = GC_M2_DAY_TO_PGC_MONTH * sum(area * value),
    .groups = 'drop'
  ) %>%
  group_by(inventory) %>%
  summarise(
    value = mean(value),
    .groups = 'drop'
  ) %>%
  mutate(estimate = 'FLUXCOM')

perturbations_annual_average <- perturbations_base %>%
  group_by(inventory, basis_vector) %>%
  summarise(
    # average over 6 years
    value = KG_M2_S_TO_PGC_MONTH * sum(area * value) / 6,
    .groups = 'drop'
  )

X_global_annual_average <- with(perturbations_annual_average, sparseMatrix(
  i = as.integer(inventory),
  j = as.integer(basis_vector),
  x = value,
  dims = c(nlevels(inventory), nlevels(basis_vector))
))

prior_emissions_annual_average <- perturbations_annual_average %>%
  group_by(inventory) %>%
  summarise(value = sum(value), .groups = 'drop') %>%
  mutate(estimate = 'Bottom-up')

posterior_emissions_annual_average_LNLGIS <- compute_posterior_mean_sd(
  prior_emissions_annual_average,
  X_global_annual_average,
  samples_LNLGIS,
  'v2.0 posterior'
)
posterior_emissions_annual_average_LNLGISSIF <- compute_posterior_mean_sd(
  prior_emissions_annual_average,
  X_global_annual_average,
  samples_LNLGISSIF,
  'v2.S posterior'
)

emissions_annual_average <- bind_rows(
  prior_emissions_annual_average,
  posterior_emissions_annual_average_LNLGIS,
  posterior_emissions_annual_average_LNLGISSIF
) %>%
  filter(inventory %in% c('bio_assim', 'bio_resp_tot')) %>%
  {
    x <- .

    bind_rows(
      x,
      x %>%
        group_by(estimate) %>%
        summarise(
          value = sum(value),
          value_samples = t(colSums(value_samples)),
          .groups = 'drop'
        ) %>%
        mutate(
          inventory = 'nee',
          value_sd = matrixStats::rowSds(value_samples)
        )
    )
  } %>%
  bind_rows(xbase_annual_average) %>%
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
  ) %>%
  select(c(inventory, estimate, value, value_sd)) %>%
  arrange(inventory, estimate)

write.csv(emissions_annual_average, args$output)
