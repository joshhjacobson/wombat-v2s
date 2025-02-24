library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(ncdf4)
library(lubridate, warn.conflicts = FALSE)
library(parallel)

source(Sys.getenv('UTILS_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--oco2-sif-directory')
parser$add_argument('--start-date')
parser$add_argument('--end-date')
parser$add_argument('--output')
args <- parser$parse_args()

sif_lite_get_date <- function(filename) {
  paste0('20', stringr::str_extract(filename, '\\d{6}'))
}

mip_average <- function(value, measurement_error) {
  sum(value * measurement_error^(-2)) / sum(measurement_error^(-2))
}

mip_error <- function(count, measurement_error) {
  1 / sqrt(sum(measurement_error^(-2)) / count)
}

log_info('Loading OCO-2 SIF observations from {args$oco2_observations_sif}')
sif_paths <- list.files(
  args$oco2_sif_directory,
  pattern = '.nc4$',
  full.names = TRUE,
  recursive = TRUE
)
sif_times <- strptime(sif_lite_get_date(basename(sif_paths)), '%Y%m%d', tz = 'UTC')
sif_paths <- sif_paths[
  sif_times >= args$start_date & sif_times < args$end_date
]
sif_soundings <- bind_rows(mclapply(sif_paths, function(filename) {
  log_trace('Loading {filename}')
  with_nc_file(list(fn = filename),
    {
      v <- function(...) as.vector(ncvar_get(fn, ...))
      oco2_soundings <- tibble(
        time = ncvar_get_time(fn, 'Delta_Time'),
        longitude = v('Longitude'),
        latitude = v('Latitude'),
        value = v('SIF_740nm'),
        measurement_error = v('SIF_Uncertainty_740nm'),
        oco2_operation_mode = v('Metadata/MeasurementMode'),
        flag = v('Quality_Flag'),
      ) %>%
        filter(
          oco2_operation_mode < 2, # keep nadir and glint only
          flag %in% c(0, 1),
          value + 3 * measurement_error > 0
        ) %>%
        select(-flag)
    },
    # HACK(jhj): ncdf4::nc_open() prints a lot of warnings since the 'SoundingId' and
    # 'FootprintId' variables have 8-byte precision in the OCO-2 SIF Lite files. We don't
    # use these variables, so we suppress the warnings.
    quiet = TRUE
  )
}, mc.cores = get_cores()))


log_info('Computing 10s averages for OCO-2 SIF observations')
time_bins_1s <- seq(
  from = as.POSIXct(args$start_date, tz = 'UTC'),
  to = as.POSIXct(args$end_date, tz = 'UTC'),
  by = '1 sec'
)
time_bins_10s <- seq(
  from = as.POSIXct(args$start_date, tz = 'UTC'),
  to = as.POSIXct(args$end_date, tz = 'UTC'),
  by = '10 sec'
)
observations_sif <- sif_soundings %>%
  # 1-second averages
  mutate(
    time = time_bins_1s[findInterval(time, time_bins_1s)]
  ) %>%
  group_by(time, oco2_operation_mode) %>%
  summarise(
    count_1s = n(),
    longitude_1s = mip_average(longitude, measurement_error),
    latitude_1s = mip_average(latitude, measurement_error),
    value_1s = mip_average(value, measurement_error),
    measurement_error_1s = mip_error(count_1s, measurement_error),
    .groups = 'drop'
  ) %>%
  # 10-second averages
  mutate(
    time = time_bins_10s[findInterval(time, time_bins_10s)] + seconds(5)
  ) %>%
  group_by(time, oco2_operation_mode) %>%
  summarise(
    count = n(),
    longitude = mip_average(longitude_1s, measurement_error_1s),
    latitude = mip_average(latitude_1s, measurement_error_1s),
    value = mip_average(value_1s, measurement_error_1s),
    measurement_error = mip_error(count, measurement_error_1s),
    .groups = 'drop'
  ) %>%
  mutate(
    observation_id = paste0(
      stringr::str_replace_all(time, "[^[:digit:]]", ""),
      oco2_operation_mode
    ),
    oco2_operation_mode = if_else(
      oco2_operation_mode == 0,
      'LN_SIF',
      'LG_SIF'
    )
  ) %>%
  select(observation_id, time, everything())

stopifnot(!any(is.na(observations_sif$value)))

log_info('Writing SIF observations to {args$output}')
fst::write_fst(observations_sif, args$output)

log_info('Done')
