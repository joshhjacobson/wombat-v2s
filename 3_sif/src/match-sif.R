library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(lubridate, warn.conflicts = FALSE)
library(ncdf4)
library(parallel)

source(Sys.getenv('UTILS_PARTIAL'))

read_inventory <- function(inventory_filename) {
  log_trace('Opening {inventory_filename}')
  inventory_fn <- nc_open(inventory_filename)

  time <- ncvar_get_time(inventory_fn, 'time')
  time_width <- diff(time)[1]
  stopifnot(all(diff(time) == time_width))
  latitude <- as.vector(ncvar_get(inventory_fn, 'lat'))
  cell_height <- get_cell_height(latitude)
  longitude <- as.vector(ncvar_get(inventory_fn, 'lon'))
  cell_width <- get_cell_width(longitude)
  sif <- ncvar_get(inventory_fn, 'sif')
  assim <- ncvar_get(inventory_fn, 'assim')

  nc_close(inventory_fn)

  list(
    time = time,
    time_width = time_width,
    latitude = latitude,
    cell_height = cell_height,
    longitude = longitude,
    cell_width = cell_width,
    sif = sif,
    assim = assim
  )
}

match_observations <- function(observations, inventory) {
  period_start <- min(inventory$time) - inventory$time_width / 2
  period_end <- max(inventory$time) + inventory$time_width / 2

  log_trace('Subsetting SIF observations to between {period_start} and {period_end}')
  observations_part <- observations %>%
    filter(
      time >= period_start,
      time < period_end
    )
  if (nrow(observations_part) == 0) {
    return(tidyr::tibble(
      observation_id = factor(),
      observation_value = numeric(),
      model_time = as.POSIXct(numeric()),
      model_longitude = numeric(),
      model_latitude = numeric(),
      model_sif = numeric(),
      model_assim = numeric()
    ))
  }

  log_trace('Computing match indices for {nrow(observations_part)} observations')
  indices <- match_grid(observations_part$longitude, observations_part$latitude, inventory)
  indices$time <- match_time(observations_part$time, inventory)

  log_trace('Subsetting SIF inventory')
  indices_array <- cbind(indices$longitude, indices$latitude, indices$time)
  sif_match <- inventory$sif[indices_array]
  assim_match <- inventory$assim[indices_array]

  observations_part %>%
    select(
      observation_id,
      observation_value = value
    ) %>%
    mutate(
      observation_id = factor(observation_id),
      model_time = inventory$time[indices$time],
      model_longitude = inventory$longitude[indices$longitude],
      model_latitude = inventory$latitude[indices$latitude],
      model_sif = sif_match,
      model_assim = assim_match
    ) %>%
    as_tibble()
}

parser <- ArgumentParser()
parser$add_argument('--oco2-observations-sif')
parser$add_argument('--inventory', nargs = '+')
parser$add_argument('--linear-models')
parser$add_argument('--output')
args <- parser$parse_args()

log_info('Loading OCO-2 SIF observations from {args$oco2_observations}')
oco2_observations <- fst::read_fst(args$oco2_observations_sif)

# SiB4 inventory terminates in 2020, so we repeat the 2020 inventory for 2021
path_final <- args$inventory[which.max(gsub('.*-(\\d+)\\.nc', '\\1', args$inventory))]
sib4_inventory_final <- read_inventory(path_final)
sib4_inventory_final$time <- sib4_inventory_final$time %m+% years(1)

control_sif <- bind_rows(
  mclapply(args$inventory, function(filename) {
    sib4_inventory <- read_inventory(filename)
    match_observations(oco2_observations, sib4_inventory)
  }, mc.cores = get_cores()),
  match_observations(oco2_observations, sib4_inventory_final)
)
stopifnot(nrow(oco2_observations) == nrow(control_sif))
stopifnot(length(setdiff(oco2_observations$observation_id, control_sif$observation_id)) == 0)

log_info('Loading fitted SIF-ASSIM models from {args$linear_models}')
linear_models <- fst::read_fst(args$linear_models)

control_sif <- control_sif %>%
  mutate(
    month = month(model_time)
  ) %>%
  inner_join(
    linear_models,
    by = c('model_longitude', 'model_latitude', 'month')
  ) %>%
  mutate(
    is_outlier = observation_value < lower_fence | observation_value > upper_fence
  ) %>%
  select(
    observation_id,
    model_time,
    model_longitude,
    model_latitude,
    value = model_sif,
    model_assim,
    intercept,
    slope,
    model_error,
    lower_fence,
    upper_fence,
    is_outlier
  )

log_info('Writing matched SIF observations to {args$output}')
fst::write_fst(control_sif, args$output)

log_info('Done')
