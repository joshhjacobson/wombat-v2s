library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(patchwork)
library(stars)

source(Sys.getenv('UTILS_PARTIAL'))
source(Sys.getenv('DISPLAY_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--model-sif-assim')
parser$add_argument('--region-sf')
parser$add_argument('--output')
args <- parser$parse_args()

model <- fst::read_fst(args$model_sif_assim)
region_sf <- readRDS(args$region_sf)
earth_bbox_sf <- rnaturalearth::ne_download(
  category = 'physical',
  type = 'wgs84_bounding_box',
  returnclass = 'sf'
)

W_to_mW  <- 1e-6
model_slope <- model %>%
  select(
    longitude = model_longitude,
    latitude = model_latitude,
    month,
    value = slope
  ) %>%
  filter(abs(latitude) != 89.5) %>%
  mutate(
    month = factor(month.abb[month], levels = month.abb),
    value = W_to_mW * value
  ) %>%
  # HACK(jhj): Remove 5 outliers in March
  filter(value < 10)


model_stars <- model_slope %>%
  group_by(longitude, latitude) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = 'drop') %>%
  arrange(longitude, latitude) %>%
  st_as_stars(dims = c('longitude', 'latitude')) %>%
  st_set_crs('WGS84') %>%
  st_transform('ESRI:54012')

model_season <- model_slope %>%
  group_by(month) %>%
  summarise(
    n = n(),
    value_q25 = quantile(value, 0.25, na.rm = TRUE),
    value_q75 = quantile(value, 0.75, na.rm = TRUE),
    value = mean(value, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    month_label = sprintf('%s\n(%d)', month, n)
  )

breaks <- seq(0, 4, by = 0.5)
limits <- c(0, 4)
base_labels <- sprintf('%.1f', breaks)
labels <- ifelse(
  abs(breaks) == max(abs(breaks)),
  sprintf(
    '%s%.1f',
    ifelse(breaks < 0, '< -', '> '),
    max(abs(breaks))
  ),
  base_labels
)

output_map <- ggplot() +
  geom_sf(data = earth_bbox_sf, fill = '#dddddd', colour = 'black') +
  geom_stars(
    data = model_stars,
    aes(fill = value)
  ) +
  geom_sf(data = region_sf, fill = NA, colour = 'grey35', linewidth = 0.1) +
  geom_segment(
    data = data.frame(y = c(-23, 23, 50)),
    mapping = aes(x = -180, y = y, xend = 180, yend = y),
    colour = 'black',
    linetype = 'dashed',
    linewidth = 0.4
  ) +
  geom_text(
    data = data.frame(
      x = c(-175, -180, -180),
      y = c(-23, 23, 50),
      label = c('23°S', '23°N', '50°N')
    ),
    mapping = aes(x = x, y = y, label = label),
    size = 8/.pt,
    nudge_y = 4
  ) +
  coord_sf(
    crs = sf::st_crs('ESRI:54012'),
    default_crs = sf::st_crs('WGS84')
  ) +
  scale_fill_binned_custom(
    'imola',
    reverse = TRUE,
    breaks = breaks,
    limits = limits,
    labels = labels,
    guide = guide_coloursteps(
      title = 'Slope',
      title.hjust = 0.5,
      label.position = 'left',
      label.theme = element_text(
        size = 7,
        colour = '#23373b',
        margin = margin(0, 0.05, 0, 0, unit = 'cm')
      ),
      barwidth = 0.5,
      barheight = 8,
      frame.colour = NA
    ),
    na.value = NA
  ) +
  labs(x = NULL, y = NULL) +
  theme(
    panel.border = element_blank(),
    legend.position = 'left',
    legend.title = element_text(
      size = 9,
      colour = '#23373b',
      margin = margin(0, 0, 0.25, 0, unit = 'cm')
    ),
    legend.margin = margin(t = 0, r = -1, b = 0, l = 0, unit = 'cm'),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'cm'),
    legend.box.spacing = unit(0, 'cm'),
    plot.title = element_blank(),
    plot.margin = margin(t = 0, b = 0, l = 0, r = 0, unit = 'cm')
  )

output_season <- model_season %>%
  ggplot(aes(x = month)) +
    geom_pointrange(
      aes(
        y = value,
        ymin = value_q25,
        ymax = value_q75
      ),
      colour = 'grey15',
      shape = 15,
      size = 0.3,
      linewidth = 0.6
    ) +
    scale_x_discrete(
      breaks = model_season$month,
      labels = model_season$month_label
    ) +
    labs(x = NULL, y = 'Slope') +
    theme(
      axis.title.y = element_text(size = 9, colour = '#23373b'),
      axis.text.y = element_text(size = 7, colour = '#23373b'),
      axis.text.x = element_text(size = 8, colour = '#23373b'),
      plot.title = element_blank(),
      plot.margin = margin(t = 0, b = 0.1, l = 0.05, r = 0.1, unit = 'cm')
    )

output <- wrap_plots(
  output_map,
  output_season,
  ncol = 1,
  heights = c(3, 1),
  widths = c(2, 1)
)

ggsave_base(
  args$output,
  output,
  width = DISPLAY_SETTINGS$full_width,
  height = 9.76
)
