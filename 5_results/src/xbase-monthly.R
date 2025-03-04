library(argparse)
library(dplyr, warn.conflicts = FALSE)

source(Sys.getenv('UTILS_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--input-files', nargs = '+')
parser$add_argument('--control-emissions')
parser$add_argument('--output')
args <- parser$parse_args()

flux_key <- c(
  'GPP' = 'bio_assim',
  'TER' = 'bio_resp_tot',
  'NEE' = 'nee'
)

read_xbase <- function(filename) {
  fn <- ncdf4::nc_open(filename)
  on.exit(ncdf4::nc_close(fn))
  v <- function(...) as.vector(ncdf4::ncvar_get(fn, ...))
  fields <- names(fn$var)
  stopifnot(sum(fields %in% names(flux_key)) == 1)
  field_name <- fields[fields %in% names(flux_key)]
  flux_name <- flux_key[field_name]
  expand.grid(
    longitude = v('lon'),
    latitude = v('lat'),
    time = ncvar_get_time(fn, 'time'),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      inventory = flux_name,
      value = v(field_name)
    )
}

cell_area <- fst::read_fst(args$control_emissions) %>%
  distinct(longitude, latitude, cell_height, area) %>%
  mutate(
    latitude_bottom = latitude - cell_height / 2
  ) %>%
  select(-cell_height)

log_debug('Reading aggregated X-Base data')
output <- lapply(args$input_files, read_xbase) %>%
  bind_rows() %>%
  mutate(
    inventory = factor(inventory, levels = c('bio_assim', 'bio_resp_tot', 'nee')),
    value = if_else(inventory == 'bio_assim', -value, value)
  ) %>%
  tidyr::drop_na() %>%
  left_join(cell_area, by = c('longitude', 'latitude'))

log_debug('Saving to {args$output}')
fst::write_fst(output, args$output)

log_debug('Done')
