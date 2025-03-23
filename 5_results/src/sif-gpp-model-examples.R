library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(ggplot2, warn.conflicts = FALSE)
library(lubridate, warn.conflicts = FALSE)
library(parallel)
library(patchwork)
library(rnaturalearth)
library(tidyr, warn.conflicts = FALSE)

source(Sys.getenv('UTILS_PARTIAL'))
source(Sys.getenv('DISPLAY_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--oco2-observations-sif')
parser$add_argument('--control-sif')
parser$add_argument('--inventory', nargs = '+')
parser$add_argument('--region-sf')
parser$add_argument('--output')
args <- parser$parse_args()


region_sf <- readRDS(args$region_sf)
earth_bbox_sf <- rnaturalearth::ne_download(
  category = 'physical',
  type = 'wgs84_bounding_box',
  returnclass = 'sf'
)
oco2_observations_sif <- fst::read_fst(args$oco2_observations_sif)
control_sif <- fst::read_fst(args$control_sif)

stopifnot(all(
  levels(oco2_observations_sif$observation_id) == levels(control_sif$observation_id)
))

control <- oco2_observations_sif %>%
  rename(value_oco2 = value) %>%
  inner_join(
    control_sif %>% rename(value_control = value),
    by = 'observation_id'
  ) %>%
  mutate(month = month(time))

inventory <- bind_rows(mclapply(args$inventory, function(path) {
  fn <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(fn))
  v <- function(...) ncdf4::ncvar_get(fn, ...)
  expand.grid(
    model_longitude = as.vector(v('lon')),
    model_latitude = as.vector(v('lat')),
    model_time = ncvar_get_time(fn, 'time'),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      sif = as.vector(v('sif')),
      assim = as.vector(v('assim')),
      month = month(model_time),
      # The below equation for local time assumes model_time is in UTC
      local_hour = hour(model_time + hours(round(model_longitude / 15)))
    ) %>%
    filter(between(local_hour, 12, 15))
}, mc.cores = get_cores())) %>% as_tibble()

control_nested <- control %>%
  nest(.by = c('model_longitude', 'model_latitude', 'month')) %>%
  rename(match_data = data) %>%
  left_join(
    inventory %>%
      nest(.by = c('model_longitude', 'model_latitude', 'month')) %>%
      rename(full_inventory = data),
    by = c('model_longitude', 'model_latitude', 'month')
  )

set.seed(20250226)
control_nested_samples <- bind_rows(
  control_nested %>%
    filter(month == 3) %>%
    sample_n(1),
  control_nested %>%
    filter(month == 7) %>%
    sample_n(1)
)

plot_scatter_location <- function(nested_data, i) {
  row <- nested_data[i, ]
  match_data <- row %>%
    select(match_data) %>%
    unnest(cols = c(match_data))
  full_inventory <- row %>%
    select(full_inventory) %>%
    unnest(cols = c(full_inventory))

  stopifnot(length(unique(match_data$intercept)) == 1)
  stopifnot(length(unique(match_data$slope)) == 1)
  intercept <- match_data$intercept[1]
  slope <- match_data$slope[1]

  p_location <- ggplot() +
    geom_sf(data = earth_bbox_sf, fill = NA, colour = 'black') +
    geom_sf(data = region_sf, fill = NA, colour = 'grey35', linewidth = 0.1) +
    coord_sf(
      crs = sf::st_crs('ESRI:54012'),
      default_crs = sf::st_crs('WGS84')
    ) +
    geom_point(
      data = row,
      aes(x = model_longitude, y = model_latitude),
      colour = '#cc3311',
      size = 2,
      shape = 15
    ) +
    labs(x = NULL, y = NULL, title = NULL) +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 0, b = 0, l = 0, r = 0, unit = 'cm')
    )

  specs <- paste0(
    'Slope: ',
    formatC(slope, format = 'e', digits = 1),
    ',  Intercept: ',
    round(intercept, 2),
    ',  Month: ',
    month.abb[row$month]
  )
  p_data <- ggplot() +
    geom_hline(
      data = match_data,
      aes(yintercept = value_oco2, colour = is_outlier),
      alpha = 0.3
    ) +
    geom_point(
      data = full_inventory,
      aes(x = assim, y = sif, shape = 'SiB4'),
      colour = '#bbbbbb',
      alpha = 0.6,
    ) +
    geom_abline(
      aes(intercept = intercept, slope = slope),
      colour = '#23373b',
      linetype = 'dashed'
    ) +
    scale_shape_manual(
      values = c('SiB4' = 1)
    ) +
    scale_colour_manual(
      values = c('FALSE' = '#009988', 'TRUE' = '#ee3377'),
      labels = c(
        'FALSE' = 'OCO-2 SIF (assimilated)',
        'TRUE' = 'OCO-2 SIF (not assimilated)'
      ),
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    expand_limits(x = 0, y = 0) +
    labs(
      x = bquote('GPP [kgCO'[2]~m^{-2}~s^{-1}*']'),
      y = bquote('SIF ['~W~m^{-2}~mu*m^{-1}~sr^{-1}*']'),
      title = specs,
      shape = NULL,
      colour = NULL
    ) +
    theme(
      plot.title = element_text(
        size = 10,
        colour = '#23373b',
        margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'cm')
      ),
      axis.title.y = element_text(size = 9, colour = '#23373b'),
      axis.text.y = element_text(size = 7, colour = '#23373b'),
      axis.title.x = element_text(size = 9, colour = '#23373b'),
      axis.text.x = element_text(size = 7, colour = '#23373b')
    )

  wrap_plots(
    p_data,
    p_location,
    nrow = 1,
    widths = c(3, 2),
    guides = 'collect'
  )
}

base_teme <- theme(
  plot.margin = margin(t = 0.1, b = 0.4, l = 0, r = 0.1, unit = 'cm'),
  legend.position = 'none'
)

bottom_theme <- theme(
  plot.margin = margin(t = 0, b = 0, l = 0, r = 0.1, unit = 'cm'),
  legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'cm'),
  legend.box.margin = margin(t = 0.2, r = 0, b = 0.1, l = -1.0, unit = 'cm'),
  legend.text = element_text(size = 9, colour = '#23373b'),
  legend.box.spacing = unit(0, 'cm'),
  legend.position = 'bottom',
  legend.box = 'horizontal'
)

p_list <- list()
for (i in seq_len(nrow(control_nested_samples))) {
  if (i == nrow(control_nested_samples)) {
    p_list[[i]] <- plot_scatter_location(control_nested_samples, i) & bottom_theme
  } else {
    p_list[[i]] <- plot_scatter_location(control_nested_samples, i) & base_teme
  }
}

output <- wrap_plots(
  p_list,
  ncol = 1,
  guides = 'collect'
)

ggsave_base(
  args$output,
  output,
  width = DISPLAY_SETTINGS$supplement_full_width,
  height = 12.5
)
