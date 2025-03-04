library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(Matrix)
library(Rcpp)
library(logger)

source(Sys.getenv('UTILS_PARTIAL'))
sourceCpp(Sys.getenv('UTILS_CPP_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--basis-vectors')
parser$add_argument('--control-emissions')
parser$add_argument('--perturbations')
parser$add_argument('--output')
args <- parser$parse_args()

basis_vectors <- fst::read_fst(args$basis_vectors)
control_emissions_bio <- fst::read_fst(args$control_emissions) %>%
  filter(inventory %in% c('bio_assim', 'bio_resp_tot'))
perturbations_bio <- fst::read_fst(args$perturbations) %>%
  filter(inventory %in% c('bio_assim', 'bio_resp_tot'))

# Construct constraint for the region total
regional_perturbations_total <- perturbations_bio %>%
  left_join(
    control_emissions_bio %>%
      distinct(longitude, latitude, area),
    by = c('longitude', 'latitude')
  ) %>%
  add_basis_vector(basis_vectors) %>%
  group_by(basis_vector, inventory, region, time) %>%
  summarise(value = sum(value * area), .groups = 'drop') %>%
  mutate(
    value = ifelse(inventory == 'bio_assim', -value, value)
  )

baseline_total <- regional_perturbations_total %>%
  group_by(inventory, region, time) %>%
  summarise(value = sum(value), .groups = 'drop') %>%
  mutate(cell_id = 1:n())

regional_perturbations_total <- regional_perturbations_total %>%
  left_join(
    baseline_total %>%
      select(inventory, region, time, cell_id),
    by = c('inventory', 'region', 'time')
  )

F_constraint_regional <- with(regional_perturbations_total, sparseMatrix(
  i = cell_id,
  j = as.integer(basis_vector),
  x = value,
  dims = c(nrow(baseline_total), nrow(basis_vectors))
))

g_constraint_regional <- baseline_total$value

# Construct constraint for the residual
with(
  basis_vectors %>%
    filter(
      inventory %in% c('bio_assim', 'bio_resp_tot'),
      component == 'residual'
    ),
  {
    baseline_residual <<- data.frame(
      region = region,
      inventory = inventory,
      longitude = NA,
      latitude = NA,
      time = NA,
      value = NA
    )
    F_residual <<- sparseMatrix(
      i = seq_along(basis_vector),
      j = as.integer(basis_vector),
      dims = c(length(basis_vector), nlevels(basis_vector))
    )
    g_residual <<- rep(1, length(basis_vector))
  }
)

log_info('Saving to {args$output}')
output <- list(
  F_sign = F_constraint_regional,
  g_sign = g_constraint_regional,
  baseline_sign = baseline_total,
  F_residual = F_residual,
  g_residual = g_residual,
  baseline_residual = baseline_residual
)
saveRDS(output, args$output)
log_info('Done')
