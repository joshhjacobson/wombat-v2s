
ar1_Q <- function(n_times, rho, sparse = TRUE) {
  stopifnot(rho >= -1 && rho <= 1)

  if (sparse) {
    if (n_times == 1) {
      return(t(sparseMatrix(i = 1, j = 1, x = 1, symmetric = TRUE)))
    }

    # Transpose ensures this is upper-triangular, the convention for this package
    t(sparseMatrix(
      i = c(
        # (1, 1) and (n_times, n_times)
        1, n_times,
        # Rest of the diagonal
        if (n_times > 2) (2 : (n_times - 1)) else NULL,
        # One off-diagonal (the other comes in via symmetry)
        2 : n_times
      ),
      j = c(
        1, n_times,
        if (n_times > 2) (2 : (n_times - 1)) else NULL,
        1 : (n_times - 1)
      ),
      x = c(
        1, 1,
        rep(1 + rho ^ 2, n_times - 2),
        rep(-rho, n_times - 1)
      ) / (1 - rho ^ 2),
      symmetric = TRUE
    ))
  } else {
    if (n_times == 1) return(matrix(1, nrow = 1))

    output <- matrix(0, nrow = n_times, ncol = n_times)
    diag(output) <- c(1, rep(1 + rho ^ 2, n_times - 2), 1)
    output[row(output) - col(output) == 1] <- -rho
    output[row(output) - col(output) == -1] <- -rho
    output / (1 - rho ^ 2)
  }
}

get_Q_block <- function(rho, w) {
  stopifnot(length(rho) == 1)
  stopifnot(length(w) == 2)
  solve(rbind(
    c(1 / w[1], rho / sqrt(prod(w))),
    c(rho / sqrt(prod(w)), 1 / w[2])
  ))
}

get_diagonal_pairs <- function(indices_list) {
  cbind(
    # NOTE(mgnb): if resp is present, the rbind cause these to alternate between assim/resp
    as.vector(do.call(rbind, indices_list)),
    as.vector(do.call(rbind, indices_list))
  )
}

get_off_diagonal_pairs <- function(indices_list) {
  rbind(
    cbind(indices_list$bio_assim, indices_list$bio_resp_tot),
    cbind(indices_list$bio_resp_tot, indices_list$bio_assim)
  )
}

get_bio_indices <- function(basis_vectors) {
  # Seasonal components
  bio_non_linear_indices <- with(basis_vectors, list(
    bio_assim = which(inventory == 'bio_assim' & !(component %in% c('intercept', 'trend', 'residual'))),
    bio_resp_tot = which(inventory == 'bio_resp_tot' & !(component %in% c('intercept', 'trend', 'residual')))
  ))
  bio_non_linear_diagonals <- get_diagonal_pairs(bio_non_linear_indices)
  bio_non_linear_off_diagonals <- get_off_diagonal_pairs(bio_non_linear_indices)

  # Residual components
  bio_residual_indices_i <- lapply(bio_regions, function(region_i) {
    which(with(basis_vectors, {
      region == region_i & inventory != 'ocean' & component == 'residual'
    }))
  })

  n_times <- length(bio_residual_indices_i[[1]]) / 2

  list(
    bio_non_linear_indices = bio_non_linear_indices,
    bio_non_linear_diagonals = bio_non_linear_diagonals,
    bio_non_linear_off_diagonals = bio_non_linear_off_diagonals,
    bio_residual_indices_i = bio_residual_indices_i,
    n_times = n_times
  )
}

get_alpha_prior_precision <- function(
  rho_bio_season,
  w_bio_season,
  kappa_bio_resid,
  rho_bio_resid,
  w_bio_resid,
  bio_indices,
  alpha_prior_precision_base
) {
  output <- alpha_prior_precision_base

  # Bio seasonal process
  Q_pair <- get_Q_block(rho_bio_season, w_bio_season)
  output[bio_indices$bio_non_linear_diagonals] <- diag(Q_pair)
  output[bio_indices$bio_non_linear_off_diagonals] <- Q_pair[1, 2]

  # Bio residual process
  Q_ar <- ar1_Q(bio_indices$n_times, kappa_bio_resid, sparse = FALSE)
  Q_cross <- get_Q_block(rho_bio_resid, w_bio_resid)
  Q_prod <- kronecker(Q_ar, Q_cross)
  for (i in seq_along(bio_regions)) {
    indices <- bio_indices$bio_residual_indices_i[[i]]
    output[indices, indices] <- Q_prod
  }

  output
}
