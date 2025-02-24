library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(sf)

source(Sys.getenv('UTILS_PARTIAL'))
source(Sys.getenv('DISPLAY_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--region-sf')
parser$add_argument('--output')
args <- parser$parse_args()

earth_bbox_sf <- rnaturalearth::ne_download(
  category = 'physical',
  type = 'wgs84_bounding_box',
  returnclass = 'sf'
)
region_sf <- readRDS(args$region_sf)

region_colours <- region_sf$colour
names(region_colours) <- region_sf$region_code

output <- ggplot() +
  geom_sf(data = earth_bbox_sf, fill = 'white', colour = 'grey35') +
  geom_sf(
    data = region_sf,
    mapping = aes(fill = region_code),
    colour = 'grey35',
    linewidth = 0.1
  ) +
  geom_sf_text(
    data = region_sf,
    mapping = aes(label = region_code),
    colour = 'black',
    nudge_x = region_sf$label_nudge_x,
    nudge_y = region_sf$label_nudge_y
  ) +
  coord_sf(
    crs = sf::st_crs('ESRI:54012'),
    default_crs = sf::st_crs('WGS84')
  ) +
  scale_fill_manual(values = region_colours) +
  guides(fill = 'none') +
  labs(x = NULL, y = NULL, title = NULL) +
  theme(
    panel.border = element_blank(),
    plot.margin = margin(t = 0, b = 0, l = 0.1, r = 0.1, unit = 'cm')
  )

ggsave_base(
  args$output,
  output,
  width = DISPLAY_SETTINGS$supplement_full_width,
  height = 8.6,
)
