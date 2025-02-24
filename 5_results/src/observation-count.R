library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(lubridate, warn.conflicts = FALSE)
library(patchwork)
library(stars)

source(Sys.getenv('UTILS_PARTIAL'))
source(Sys.getenv('DISPLAY_PARTIAL'))


match_grid <- function(df, model) {
  latitude_index <- findInterval(
    df$latitude,
    c(model$latitude[1] - model$cell_height[1] / 2, model$latitude + model$cell_height / 2),
    rightmost.closed = TRUE
  )
  longitude_index <- findInterval(
    df$longitude,
    c(model$longitude[1] - model$cell_width[1] / 2, model$longitude + model$cell_width / 2),
    rightmost.closed = FALSE
  )
  longitude_index[longitude_index > length(model$longitude)] <- 1L

  df %>%
    mutate(
      latitude = model$latitude[latitude_index],
      longitude = model$longitude[longitude_index]
    )
}


parser <- ArgumentParser()
parser$add_argument('--observations')
parser$add_argument('--control', nargs = '+')
parser$add_argument('--region-sf')
parser$add_argument('--output')
args <- parser$parse_args()


earth_bbox_sf <- rnaturalearth::ne_download(
  category = 'physical',
  type = 'wgs84_bounding_box',
  returnclass = 'sf'
)
region_sf <- readRDS(args$region_sf)

control <- bind_rows(lapply(args$control, fst::read_fst)) %>%
  mutate(observation_id = droplevels(observation_id))

# Get the assimilated observations
observations <- fst::read_fst(args$observations) %>%
  filter(
    assimilate %in% c(1, 3),
    observation_id %in% control$observation_id,
    overall_observation_mode != 'OG'
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

observations <- bind_rows(
  observations %>%
    filter(observation_type != 'oco2'),
  observations %>%
    filter(observation_type == 'oco2') %>%
    mutate(
      observation_type = if_else(
        overall_observation_mode %in% c('LN_SIF', 'LG_SIF'),
        'oco2_sif',
        'oco2_xco2'
      )
    )
)


grid_system <- list(
  latitude = seq(-88, 88, by = 2),
  cell_height = 2,
  longitude = seq(-180, 177.5, by = 2.5),
  cell_width = 2.5
)
grid_base <- expand.grid(
  latitude = grid_system$latitude,
  longitude = grid_system$longitude
)

observation_types <- c('oco2_sif', 'oco2_xco2', 'obspack')
strip_labels <- expression(
  oco2_sif = 'OCO-2 SIF',
  oco2_xco2 = 'OCO-2 XCO'[2],
  obspack = 'In-situ/flask mole fraction'
)

output_columns <- lapply(seq_along(observation_types), function(i) {
  count_space <- grid_base %>%
    left_join(
      observations %>%
        filter(
          observation_type == observation_types[i],
          abs(latitude) <= 89
        ) %>%
        match_grid(grid_system) %>%
        group_by(longitude, latitude) %>%
        summarise(n = n(), .groups = 'drop'),
      by = c('longitude', 'latitude')
    ) %>%
    arrange(longitude, latitude) %>%
    st_as_stars(dims = c('longitude', 'latitude')) %>%
    st_set_crs('WGS84') %>%
    st_transform('ESRI:54012')

  count_time <- observations %>%
    filter(observation_type == observation_types[i]) %>%
    mutate(
      year_month = as.Date(floor_date(time, 'month') + days(15))
    ) %>%
    group_by(year_month) %>%
    summarise(n = 1e-3 * n(), .groups = 'drop')

  output_space <- ggplot() +
    geom_sf(data = earth_bbox_sf, fill = 'grey85', colour = 'black') +
    geom_stars(
      data = count_space,
      aes(fill = n)
    ) +
    geom_sf(data = region_sf, fill = NA, colour = 'grey35', linewidth = 0.1) +
    coord_sf(
      crs = sf::st_crs('ESRI:54012'),
      default_crs = sf::st_crs('WGS84')
    ) +
    scale_fill_binned_custom(
      'glasgow',
      reverse = TRUE,
      transform = 'log10',
      breaks = 10 ^ seq(0, 5, by = 0.5),
      limits = c(1, 10^5),
      labels = function(x) {
        ifelse((x %% 1 == 0) & (x > 1) & (x < 10^5), sprintf('%.0f', x), '')
      },
      guide = guide_coloursteps(
        title.position = 'left',
        axis = FALSE,
        label.theme = element_text(
          size = 7,
          colour = '#23373b',
          margin = margin(0.1, 0, 0.1, 0, unit = 'cm')
        ),
        frame.colour = NA,
        barwidth = 13,
        barheight = 0.5
      ),
      na.value = NA
    ) +
    labs(x = NULL, y = NULL, title = NULL, fill = 'Obs. count') +
    theme(
      panel.border = element_blank()
    )

  output_time <- ggplot(count_time, aes(x = year_month, y = n)) +
    geom_line(colour = 'grey15') +
    scale_x_date(date_labels = '%Y-%m') +
    labs(
      x = NULL,
      y = NULL,
      title = strip_labels[i]
    ) +
    theme(
      plot.title = element_text(size = 9, hjust = 0.5, margin = margin(0.1, 0, 0, 0, unit = 'cm')),
      axis.title.y = element_text(size = 8),
      axis.text.y = element_text(size = 7),
      axis.text.x = element_text(size = 7),
      axis.title.x = element_blank()
    )
  if (i == 1) {
    output_time <- output_time + labs(y = 'Obs. count [K]')
  }

  wrap_plots(
    output_time,
    output_space,
    ncol = 1,
    heights = c(1, 2),
    widths = c(1, 1)
  )
})

output <- wrap_plots(output_columns, ncol = 3, guides = 'collect') & theme(
  plot.margin = margin(t = 0, r = 0.1, b = 0, l = 0.05, unit = 'cm'),
  legend.position = 'bottom',
  legend.title = element_text(size = 8, colour = '#23373b', vjust = 1),
  legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'cm'),
  legend.box.margin = margin(t = -0.25, r = 0, b = 0, l = -1.5, unit = 'cm')
)

ggsave_base(
  args$output,
  output,
  width = DISPLAY_SETTINGS$supplement_full_width,
  height = 5.25
)
