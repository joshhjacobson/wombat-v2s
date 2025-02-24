library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(Matrix)
library(Rcpp)

rcpp_cache_dir <- Sys.getenv('RCPP_CACHE_DIR')
options(rcpp.cache.dir = if (rcpp_cache_dir == '') tempdir() else rcpp_cache_dir)

source(Sys.getenv('UTILS_PARTIAL'))
source(Sys.getenv('ALPHA_PRECISION_PARTIAL'))
sourceCpp(Sys.getenv('UTILS_CPP_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--prior-base')
parser$add_argument('--basis-vectors')
parser$add_argument('--control-emissions')
parser$add_argument('--perturbations')
parser$add_argument('--output')
args <- parser$parse_args()


prior_base <- readRDS(args$prior_base)
basis_vectors <- fst::read_fst(args$basis_vectors)
control_emissions <- fst::read_fst(args$control_emissions)
perturbations_base <- fst::read_fst(args$perturbations)

cell_area <- control_emissions %>%
  distinct(longitude, latitude, cell_height, area) %>%
  mutate(
    latitude_bottom = latitude - cell_height / 2
  ) %>%
  select(-cell_height)

perturbations <- perturbations_base %>%
  add_basis_vector(basis_vectors) %>%
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
  ) %>%
  left_join(cell_area, by = c('longitude', 'latitude'))

# NOTE(mgnb): values of alpha fixed to zero are given infinite precision
alpha_to_include <- is.finite(diag(prior_base$precision))

alpha_prior_precision_base <- prior_base$precision[alpha_to_include, alpha_to_include]
n_alpha <- nrow(alpha_prior_precision_base)

bio_regions <- sort(unique(basis_vectors$region))[1 : 11]
bio_indices <- get_bio_indices(basis_vectors[alpha_to_include, ])

get_emissions_parts <- function(filter_expr) {
  perturbations_region <- perturbations %>%
    filter({{ filter_expr }}) %>%
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
    mutate(output = 'Bottom-up')

  list(
    X_region = X_region,
    prior_emissions = prior_emissions
  )
}

get_linear_transform <- function(region_name, inventory_name) {
  emission_parts <- get_emissions_parts(
    region == region_name & inventory == inventory_name & minor_component == 'linear'
  )
  linear_columns <- which(with(
    basis_vectors[alpha_to_include, ],
    region == region_name & inventory == inventory_name & component %in% c('intercept', 'trend')
  ))

  X_region_linear <- emission_parts$X_region[, linear_columns]
  intercept <- X_region_linear[1, 1]
  mean_trend <- mean(X_region_linear[, 2])
  R_block <- rbind(
    c(1 + mean_trend / intercept, -mean_trend / intercept),
    c(0, 1)
  )

  list(
    indices = linear_columns,
    R_block = R_block
  )
}

get_transformation_matrix <- function(region_name) {
  transform_bio_assim <- get_linear_transform(region_name, 'bio_assim')
  transform_bio_resp_tot <- get_linear_transform(region_name, 'bio_resp_tot')

  R_transform <- Matrix::bdiag(
    transform_bio_assim$R_block,
    transform_bio_resp_tot$R_block
  )

  list(
    R_transform = R_transform,
    indices = c(transform_bio_assim$indices, transform_bio_resp_tot$indices)
  )
}

log_debug('Computing transformation matrix')
R_transform <- Diagonal(n_alpha)
for (region_name in sprintf('Region%02d', 1 : 11)) {
  region_transform <- get_transformation_matrix(region_name)
  R_transform[
    region_transform$indices,
    region_transform$indices
  ] <- region_transform$R_transform
}

Q_base <- get_alpha_prior_precision(
  0,  # rho season
  c(1, 1),  # w season
  0,  # kappa resid
  0,  # rho resid
  c(1, 1),  # w resid
  bio_indices,
  alpha_prior_precision_base
)

trend_indices <- which(with(basis_vectors[alpha_to_include, ], component == 'trend'))
Q_base[cbind(trend_indices, trend_indices)] <- 1e-4
Q <- as.matrix(solve(t(R_transform), t(solve(t(R_transform), Q_base))))

alpha_prior_precision_final <- prior_base$precision
alpha_prior_precision_final[alpha_to_include, alpha_to_include] <- Q

log_debug('Saving to {args$output}')
saveRDS(list(
  mean = prior_base$mean,
  precision = alpha_prior_precision_final
), args$output)

log_debug('Done')
