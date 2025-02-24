library(logger)

if (Sys.getenv('WOMBAT_LOG_LEVEL') != '') {
  log_threshold(list(
    'trace' = TRACE,
    'debug' = DEBUG,
    'info' = INFO,
    'warn' = WARN,
    'error' = ERROR
  )[[Sys.getenv('WOMBAT_LOG_LEVEL')]])
}

get_cores <- function() {
  omp_threads <- as.integer(Sys.getenv('OMP_NUM_THREADS'))
  max_workers <- as.integer(Sys.getenv('WOMBAT_MAX_WORKERS'))
  if (!is.na(omp_threads) & omp_threads > 0) {
    omp_threads
  } else if (!is.na(max_workers) & max_workers > 0) {
    max_workers
  } else {
    parallel::detectCores()
  }
}

# HACK(mgnb): as.integer fails when any factor entry is NA; this preserves
# the NA values
to_integer <- function(x) {
  class(x) <- NULL
  attr(x, 'levels') <- NULL
  x
}

add_basis_vector <- function(x, basis_vectors) {
  for (name in c('inventory', 'component', 'region', 'month')) {
    if (!identical(levels(x[[name]]), levels(basis_vectors[[name]]))) {
      x[[name]] <- update_levels(x[[name]], levels(basis_vectors[[name]]))
    }
  }

  x %>%
    mutate(
      basis_vector = to_basis_vector(
        to_integer(inventory),
        to_integer(component),
        to_integer(region),
        to_integer(month),
        to_integer(basis_vectors$inventory),
        to_integer(basis_vectors$component),
        to_integer(basis_vectors$region),
        to_integer(basis_vectors$month),
        levels(basis_vectors$basis_vector)
      )
    )
}

process_sensitivities <- function(
  control_path,
  sensitivities_base,
  basis_vectors,
  observation_id_levels,
  fn
) {
  load_part <- function(path) {
    log_debug('Reading {path}')
    output <- fst::read_fst(path)

    log_trace('Filtering')
    if (grepl('daily', path)) {
      output <- output %>%
        filter(
          is_level_in(observation_id, observation_id_levels),
          component == 'residual'
        )
    } else {
      output <- output %>%
        filter(
          is_level_in(observation_id, observation_id_levels)
        )
    }

    log_trace('Releveling')
    output <- output %>%
      mutate(
        observation_id = update_levels(
          observation_id,
          observation_id_levels
        )
      )

    log_trace('Adding basis_vector')
    output %>%
      add_basis_vector(basis_vectors) %>%
      select(observation_id, basis_vector, value)
  }

  hourly_residual_complete <- NULL
  hourly_base <- stringr::str_replace(basename(control_path), '\\.fst', '')
  hourly_paths <- list.files(
    dirname(sensitivities_base),
    pattern = sprintf('%s-%s.*', basename(sensitivities_base), hourly_base),
    full.names = TRUE
  )
  for (path in hourly_paths) {
    sensitivities <- load_part(path)

    fn(sensitivities, rep(FALSE, nrow(sensitivities)))

    log_trace('Archiving residual parts')
    hourly_residual_complete <- c(
      hourly_residual_complete,
      list(
        sensitivities %>%
          select(observation_id, basis_vector) %>%
          filter(basis_vector %in% basis_vectors$basis_vector[basis_vectors$component == 'residual']) %>%
          arrange(observation_id, basis_vector)
      )
    )

    log_trace('Cleaning up')
    rm(sensitivities)
    gc()
  }

  daily_base <- stringr::str_replace(hourly_base, 'hourly', 'daily')
  daily_paths <- list.files(
    dirname(sensitivities_base),
    pattern = sprintf('%s-%s.*', basename(sensitivities_base), daily_base),
    full.names = TRUE
  )
  for (path in daily_paths) {
    sensitivities <- load_part(path)

    log_trace('Finding pre-existing parts')
    is_done <- rep(FALSE, nrow(sensitivities))
    for (part in hourly_residual_complete) {
      is_done <- is_done | sensitivities_already_done(
        sensitivities$observation_id,
        sensitivities$basis_vector,
        part$observation_id,
        part$basis_vector
      )
    }

    fn(sensitivities, is_done)

    log_trace('Cleaning up')
    rm(sensitivities)
    gc()
  }

  invisible(NULL)
}

ou_precision <- function(diff_x, lambda, precisions, cross_precisions) {
  if (length(precisions) == 1) {
    return(TridiagonalMatrix(
      precisions,
      numeric(0)
    ))
  }

  stopifnot(all(diff_x > 0))

  e_lambda_d_x <- exp(-lambda * diff_x)
  e_2lambda_d_x <- exp(-2 * lambda * diff_x)

  r <- e_2lambda_d_x / (1 - e_2lambda_d_x)
  major <- c(r, 0) + c(0, r) + 1
  minor <- -e_lambda_d_x / (1 - e_2lambda_d_x)

  TridiagonalMatrix(
    major * precisions,
    minor * cross_precisions
  )
}

chol_solve <- function(R, b) {
  backsolve(R, backsolve(R, b, transpose = TRUE))
}

ncvar_get_time <- function(x, variable_name) {
  units <- x$var[[variable_name]]$units
  if (is.null(units)) {
    units <- x$dim[[variable_name]]$units
  }

  match <- stringr::str_match(
    units,
    '^(days|hours|minutes|seconds|Days|Hours|Seconds) since (.+)$'
  )[1, 2 : 3]
  lubridate::parse_date_time(
    match[2],
    orders = c('ymd', 'ymd HM', 'ymd HMS')
  ) + lubridate::period(sprintf(
    '%.04f %s',
    ncdf4::ncvar_get(x, variable_name),
    tolower(match[1])
  ))
}

with_nc_file <- function(files, code, envir = parent.frame(), quiet = FALSE) {
  for (nme in names(files)) {
    ifelse(
      quiet,
      invisible(capture.output(assign(
        nme,
        ncdf4::nc_open(files[[nme]]),
        envir = envir
      ))),
      assign(
        nme,
        ncdf4::nc_open(files[[nme]]),
        envir = envir
      )
    )
  }
  on.exit({
    for (nme in names(files)) {
      ncdf4::nc_close(get(nme, envir = envir))
      rm(list = nme, envir = envir)
    }
  })
  eval(substitute(code), envir = envir)
}

get_cell_width <- function(longitude) {
  rep(longitude[2] - longitude[1], length(longitude))
}

get_cell_height <- function(latitude) {
  unique_diff <- unique(diff(latitude))
  if (length(unique_diff) == 1) {
    rep(unique_diff, length(latitude))
  } else {
    # Half-height polar cells
    full_width <- max(unique_diff)
    c(full_width / 2, rep(full_width, length(latitude) - 2), full_width / 2)
  }
}

match_time <- function(time, model) {
  bins <- c(model$time[1] - model$time_width / 2, model$time + model$time_width / 2)
  if (is.unsorted(bins)) {
    bins <- sort(bins)
  }
  findInterval(time, bins, rightmost.closed = FALSE)
}

match_grid <- function(longitude, latitude, model) {
  latitude_index <- findInterval(
    latitude,
    c(model$latitude[1] - model$cell_height[1] / 2, model$latitude + model$cell_height / 2),
    rightmost.closed = TRUE
  )
  longitude_index <- findInterval(
    longitude,
    c(model$longitude[1] - model$cell_width[1] / 2, model$longitude + model$cell_width / 2),
    rightmost.closed = TRUE
  )
  longitude_index[longitude_index > length(model$longitude)] <- 1L

  data.frame(
    latitude = latitude_index,
    longitude = longitude_index
  )
}

read_gridded_data <- function(
  filename,
  field_name,
  include_time = TRUE,
  include_variable = FALSE
) {
  with_nc_file(list(fn = filename), {
    output <- list(
      latitude = as.vector(ncdf4::ncvar_get(fn, 'lat')),
      longitude = as.vector(ncdf4::ncvar_get(fn, 'lon')),
      value = ncdf4::ncvar_get(fn, field_name),
      value_units = fn$var[[field_name]]$units,
      value_longname = fn$var[[field_name]]$longname
    )
    if (include_time) {
      output$time <- ncvar_get_time(fn, 'time')
    }
    if (include_variable) {
      output$variable <- as.vector(ncdf4::ncvar_get(fn, 'variable'))
    }
    attrs <- ncdf4::ncatt_get(fn, 0)
    for (name in names(attrs)) {
      attr(output, name) <- attrs[[name]]
    }
  })
  output
}

climatology_design_matrix <- function(
  time,
  origin,
  time_resolution = c('hourly', 'monthly'),
  time_alignment = c('start', 'middle', 'end'),
  harmonics
) {
  time_resolution <- match.arg(time_resolution)
  time_alignment <- match.arg(time_alignment)

  if (time_alignment != 'middle') {
    if (time_resolution == 'monthly') {
      time <- time + seconds(
        c('start' = 1, 'end' = -1)[time_alignment]
        * 24 * 60 * 60 * days_in_month(time) / 2
      )
    } else {
      time <- time + seconds(
        c('start' = 1, 'end' = -1)[time_alignment]
        * 60 * 30
      )
    }
  }

  trend <- as.numeric(
    time - origin,
    unit = 'hours'
  ) / (365.25 * 24)
  X <- matrix(NA, nrow = length(time), ncol = 4 * harmonics + 2)
  X[, 1] <- 1
  X[, 2] <- trend
  for (k in seq_len(harmonics)) {
    X[, 3 + 2 * (k - 1)] <- cospi(2 * trend / (1 / k))
    X[, 3 + 2 * (k - 1) + 1] <- sinpi(2 * trend / (1 / k))
    X[, 3 + 2 * harmonics + 2 * (k - 1)] <- trend * cospi(
      2 * trend / (1 / k)
    )
    X[, 3 + 2 * harmonics + 2 * (k - 1) + 1] <- trend * sinpi(
      2 * trend / (1 / k)
    )
  }
  colnames(X) <- c(
    'intercept',
    'trend',
    as.vector(sapply(seq_len(harmonics), function(k) {
      sprintf(c('cos12_%d_intercept', 'sin12_%d_intercept'), k)
    })),
    as.vector(sapply(seq_len(harmonics), function(k) {
      sprintf(c('cos12_%d_trend', 'sin12_%d_trend'), k)
    }))
  )
  X
}

read_climatology <- function(filename) {
  gridded <- read_gridded_data(
    filename,
    'coefficient',
    include_time = FALSE,
    include_variable = TRUE
  )
  expand.grid(
    longitude = gridded$longitude,
    latitude = gridded$latitude,
    variable = gridded$variable,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      variable = factor(variable, gridded$variable),
      value = as.double(gridded$value)
    )
}

discretise_by_breaks <- function(y, breaks, limits) {
  full_breaks <- c(limits[1], breaks, limits[2])
  midpoints <- (head(full_breaks, -1) + tail(full_breaks, -1)) / 2
  midpoints[cut(y, c(-Inf, breaks, Inf))]
}

compute_posterior <- function(prior, design_matrix, samples, posterior_name) {
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
      value_q025 = matrixStats::rowQuantiles(value_samples, probs = 0.025),
      value_q975 = matrixStats::rowQuantiles(value_samples, probs = 0.975)
    ) %>%
    select(-value_prior)
}

# These are from the WOMBATv1 paper
W_PRIOR <- list(
  shape = 0.3542832,
  rate = 0.01534294
)

GAMMA_PRIOR <- list(
  shape = 1.62702,
  rate = 2.171239
)

ELL_PRIOR <- list(
  shape = 1,
  rate = 1
)

N_MCMC_SAMPLES <- 5000
N_MCMC_WARM_UP <- 1000

PER_SECONDS_TO_PER_YEAR <- 365.25 * 24 * 60 * 60
KG_M2_S_TO_PGC_MONTH <- 12.01 / 44.01 * 1e-12 * 365.25 * 24 * 60 * 60 / 12
KG_M2_S_TO_PGC_YEAR <- 12.01 / 44.01 * 1e-12 * 365.25 * 24 * 60 * 60
GC_M2_DAY_TO_PGC_MONTH <- 1e-15 * 365.25 / 12
GC_M2_DAY_TO_PGC_YEAR <- 1e-15 * 365.25
GC_DAY_TO_KGCO2_YEAR <- 44.01 / 12.01 * 1e-3 * 365.25
