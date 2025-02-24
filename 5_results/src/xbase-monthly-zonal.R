library(argparse)
library(dplyr, warn.conflicts = FALSE)

source(Sys.getenv('UTILS_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--xbase-monthly-2x25')
parser$add_argument('--area-1x1')
parser$add_argument('--output')
args <- parser$parse_args()

xbase_monthly_2x25 <- fst::read_fst(args$xbase_monthly_2x25)

with_nc_file(list(fn = args$area_1x1), {
  v <- function(...) as.vector(ncdf4::ncvar_get(fn, ...))
  area_1x1 <- expand.grid(
    longitude = v('lon'),
    latitude = v('lat')
  ) %>%
    mutate(area = v('cell_area'))
})

area_495 <- (area_1x1 %>% filter(latitude == 49.5) %>% pull(area))[1]
area_505 <- (area_1x1 %>% filter(latitude == 50.5) %>% pull(area))[1]
area_both <- area_495 + area_505

# Splits grid cells that cross boundaries
xbase_monthly_2x25_zonal <- bind_rows(
  xbase_monthly_2x25 %>%
    filter(latitude != 0 & latitude != 50),
  xbase_monthly_2x25 %>%
    filter(latitude == 0) %>%
    mutate(
      latitude = -0.5,
      latitude_bottom = -1,
      area = area / 2
    ),
  xbase_monthly_2x25 %>%
    filter(latitude == 0) %>%
    mutate(
      latitude = 0.5,
      latitude_bottom = 0,
      area = area / 2
    ),
  xbase_monthly_2x25 %>%
    filter(latitude == 50) %>%
    mutate(
      latitude = 49.5,
      latitude_bottom = 49,
      area = area * area_495 / area_both
    ),
  xbase_monthly_2x25 %>%
    filter(latitude == 50) %>%
    mutate(
      latitude = 50.5,
      latitude_bottom = 50,
      area = area * area_505 / area_both
    )
)

fst::write_fst(xbase_monthly_2x25_zonal, args$output)
