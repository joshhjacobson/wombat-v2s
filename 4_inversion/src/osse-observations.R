library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(fastsparse)
library(Matrix)
library(fst)

source(Sys.getenv('UTILS_PARTIAL'))

sample_normal_precision <- function(Q) {
  as.vector(solve(chol(Q), rnorm(ncol(Q))))
}

parser <- ArgumentParser()
parser$add_argument('--seed', type = 'integer', default = 0)
parser$add_argument('--true-alpha')
parser$add_argument('--basis-vectors')
parser$add_argument('--hyperparameter-estimates')
parser$add_argument('--observations')
parser$add_argument('--overall-observation-mode', nargs = '+')
parser$add_argument('--control', nargs = '+')
parser$add_argument('--component-name', nargs = '+')
parser$add_argument('--component-parts', nargs = '+')
parser$add_argument('--component-transport-matrix', nargs = '+')
parser$add_argument('--output')
args <- parser$parse_known_args()[[1]]

set.seed(20240403 + args$seed)

log_debug('Loading control')
control <- bind_rows(lapply(args$control, read_fst)) %>%
  mutate(observation_id = droplevels(observation_id))

log_debug('Loading observations')
observations <- read_fst(args$observations) %>%
  filter(
    assimilate %in% c(1, 3),
    overall_observation_mode %in% args$overall_observation_mode,
    observation_id %in% control$observation_id
  ) %>%
  arrange(observation_group, time) %>%
  left_join(
    control %>%
      select(observation_id, value_control = value),
    by = 'observation_id'
  ) %>%
  mutate(
    observation_id = factor(
      as.character(observation_id),
      as.character(observation_id)
    )
  )

control <- control %>%
  filter(observation_id %in% observations$observation_id) %>%
  mutate(
    observation_id = factor(
      as.character(observation_id),
      levels(observations$observation_id)
    )
  )

stopifnot(all(
  levels(observations$observation_id) == levels(control$observation_id)
))
stopifnot(nlevels(observations$observation_id) == nrow(observations))
stopifnot(!anyNA(observations$error))

part_indices <- seq_along(args$component_name)
observations$component_name <- ''
for (part_i in part_indices) {
  observations$component_name <- ifelse(
    observations$overall_observation_mode %in% strsplit(
      args$component_parts[part_i],
      '|',
      fixed = TRUE
    )[[1]],
    args$component_name[part_i],
    observations$component_name
  )
}

hyperparameter_groups <- observations %>%
  group_by(hyperparameter_group) %>%
  group_modify(~ tibble(
    component_names = list(unique(.x$component_name))
  ))
hyperparameter_group_indices <- seq_len(nrow(hyperparameter_groups))

observation_parts <- lapply(hyperparameter_group_indices, function(i) {
  component_names_i <- hyperparameter_groups$component_names[[i]]
  hyperparameter_group_i <- hyperparameter_groups$hyperparameter_group[i]

  lapply(component_names_i, function(component_name_i_j) {
    observations %>%
      filter(
        component_name == component_name_i_j,
        hyperparameter_group == hyperparameter_group_i
      )
  }) %>%
    bind_rows()
})

hyperparameter_estimates <- read_fst(args$hyperparameter_estimates)

log_debug('Simulating observation errors')
epsilon_parts <- lapply(hyperparameter_group_indices, function(i) {
  observations_i <- observation_parts[[i]]

  epsilon_parts_i <- observations_i %>%
    group_by(observation_group, hyperparameter_group) %>%
    group_map(~ {
      parameters <- .y %>%
        left_join(hyperparameter_estimates, by = 'hyperparameter_group')

      diff_time <- diff(as.double(
        .x$time - .x$time[1],
        unit = parameters$ell_unit
      ))
      # HACK(mgnb): some IS sites have repeated times; separate them minimally
      diff_time[diff_time == 0] <- 1 / 24

      n <- nrow(.x)
      output <- if (parameters$rho == 1) {
        sample_normal_precision(ou_precision(
          diff_time,
          1 / parameters$ell,
          rep(1, n),
          rep(1, n - 1)
        ))
      } else {
        (
          sample_normal_precision(Diagonal(
            x = rep(1 / (1 - parameters$rho), n)
          ))
          + sample_normal_precision(ou_precision(
            diff_time,
            1 / parameters$ell,
            rep(1 / parameters$rho, n),
            rep(1 / parameters$rho, n - 1)
          ))
        )
      }

      output * .x$error / sqrt(parameters$gamma)
    })

  epsilons_i <- do.call(c, epsilon_parts_i)

  stopifnot(nrow(observations_i) == length(epsilons_i))
  observations_i %>%
    select(observation_id) %>%
    mutate(
      epsilon = epsilons_i
    )
})

osse_epsilon <- bind_rows(epsilon_parts)
stopifnot(!anyNA(osse_epsilon$epsilon))
stopifnot(nrow(observations) == nrow(osse_epsilon))
stopifnot(all(
  levels(observations$observation_id) == levels(osse_epsilon$observation_id)
))

osse_observations <- observations %>%
  left_join(
    osse_epsilon,
    by = 'observation_id'
  ) %>%
  mutate(
    value = value_control + epsilon
  )

if (!is.null(args$true_alpha)) {
  log_debug('Computing perturbations with true alpha from {args$true_alpha}')
  true_alpha <- read_fst(args$true_alpha)
  basis_vectors <- read_fst(args$basis_vectors)

  n_all_alpha <- nrow(basis_vectors)

  log_trace('Loading transport matrices')
  H_parts <- lapply(part_indices, function(i) {
    log_trace('Loading {args$component_transport_matrix[i]}')
    observations_i <- observations %>%
      filter(component_name == args$component_name[i])
    n_observations_i <- observations_i %>%
      nrow() %>%
      as.double()

    n_i <- n_observations_i * n_all_alpha
    log_trace('Total size to load is {n_i} = {n_observations_i} x {n_all_alpha}')

    fn <- pipe(sprintf('lz4 -v %s -', args$component_transport_matrix[i]), 'rb')
    H_vec <- readBin(fn, 'double', n_i)
    close(fn)
    output <- matrix(H_vec, nrow = n_observations_i)
    gc()
    output[, as.integer(true_alpha$basis_vector)]
  })
  gc()

  log_trace('Adding perturbations to observations')
  osse_observations <- bind_rows(lapply(part_indices, function(i) {
    osse_observations %>%
      filter(
        component_name == args$component_name[i]
      ) %>%
      mutate(
        value = value + as.vector(H_parts[[i]] %*% true_alpha$value)
      )
  }))
} else {
  log_debug('No perturbation when true alpha is zero')
}

osse_observations <- osse_observations %>%
  select(-c(
    component_name,
    value_control
  ))

log_debug('Saving to {args$output}')
write_fst(osse_observations, args$output)

log_debug('Done')
