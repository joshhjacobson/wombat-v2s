library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(Matrix)
library(patchwork)
library(stars)

source(Sys.getenv('UTILS_PARTIAL'))
source(Sys.getenv('DISPLAY_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--six-year-average')
parser$add_argument('--region-sf')
parser$add_argument('--output')
args <- parser$parse_args()

flux_key <- list(
  name = 'bio_assim',
  label_precision = 0,
  drop_second_labels = FALSE,
  palette = 'bamako',
  reverse = FALSE,
  symmetric = FALSE,
  mean_breaks = round(seq(-12, 0, by = 2), 1),
  mean_limits = c(-14, 0),
  diff_breaks = round(seq(-6, 6, by = 2), 1),
  diff_limits = c(-8, 8),
  sd_breaks = round(seq(0, 0.4, by = 0.1), 1),
  sd_limits = c(0, 0.4 + 1.5 * 0.1)
)

region_sf <- readRDS(args$region_sf)
earth_bbox_sf <- rnaturalearth::ne_download(
  category = 'physical',
  type = 'wgs84_bounding_box',
  returnclass = 'sf'
)

six_year_average_stars <- fst::read_fst(args$six_year_average) %>%
  filter(
    inventory == flux_key$name,
    abs(latitude) != 89.5
  ) %>%
  select(-inventory) %>%
  arrange(longitude, latitude, estimate) %>%
  st_as_stars(dims = c('longitude', 'latitude', 'estimate'))

flux_label <- bquote('GPP flux [kgCO'[2]~m^{-2}~yr^{-1}*']')
diff_label <- bquote('GPP difference [kgCO'[2]~m^{-2}~yr^{-1}*']')

average_posterior_mean_wombat <- six_year_average_stars %>%
  filter(estimate == 'LNLGISSIF') %>%
  st_set_crs('WGS84') %>%
  st_transform('ESRI:54012') %>%
  plot_map(
    value,
    flux_key$mean_breaks,
    flux_key$mean_limits,
    flux_key$palette,
    reverse = flux_key$reverse,
    symmetric = flux_key$symmetric,
    label_precision = flux_key$label_precision,
    drop_second_labels = FALSE,
    show_excess = TRUE,
    bar_width = 14
  ) +
    labs(fill = flux_label, title = 'WOMBAT v2.S post. mean')

average_posterior_sd_wombat <- six_year_average_stars %>%
  filter(estimate == 'LNLGISSIF') %>%
  st_set_crs('WGS84') %>%
  st_transform('ESRI:54012') %>%
  plot_map(
    value_scale,
    flux_key$sd_breaks,
    flux_key$sd_limits,
    'devon',
    symmetric = FALSE,
    reverse = TRUE,
    trim_colours = TRUE,
    drop_second_labels = FALSE,
    label_precision = ifelse(flux_key$name == 'nee', 3, 1),
    show_excess = TRUE
  ) +
    labs(fill = flux_label, title = 'WOMBAT v2.S post. std. dev.')

average_fluxcom <- six_year_average_stars %>%
  filter(estimate == 'FLUXCOM') %>%
  st_set_crs('WGS84') %>%
  st_transform('ESRI:54012') %>%
  plot_map(
    value,
    flux_key$mean_breaks,
    flux_key$mean_limits,
    flux_key$palette,
    reverse = flux_key$reverse,
    symmetric = flux_key$symmetric,
    label_precision = flux_key$label_precision,
    drop_second_labels = FALSE,
    show_excess = TRUE,
    bar_width = 14
  ) +
    labs(fill = flux_label, title = 'FLUXCOM')

average_diff_fluxcom <- six_year_average_stars %>%
  filter(estimate == 'FLUXCOM Difference') %>%
  st_set_crs('WGS84') %>%
  st_transform('ESRI:54012') %>%
  plot_map(
    value,
    flux_key$diff_breaks,
    flux_key$diff_limits,
    'vik',
    symmetric = TRUE,
    label_precision = min(c(1, flux_key$label_precision)),
    drop_second_labels = flux_key$drop_second_labels,
    show_excess = FALSE
  ) +
    labs(fill = diff_label, title = 'WOMBAT v2.S - FLUXCOM')

base_theme <- theme(
  plot.title = element_text(
    hjust = 0.5,
    size = 9,
    margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'cm')
  ),
  plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'cm')
)

top_theme <- base_theme +
  theme(
    plot.margin = margin(t = 0, r = 0, b = 0.2, l = 0, unit = 'cm')
  )

output <- wrap_plots(
  wrap_plots(
    average_posterior_mean_wombat + top_theme,
    average_fluxcom + top_theme,
    ncol = 2,
    guides = 'collect'
  ) &
    theme(
      legend.position = 'bottom',
      legend.margin = margin(t = -0.4, r = 0, b = -0.2, l = 0, unit = 'cm')
    ),
  wrap_plots(
    average_posterior_sd_wombat + base_theme,
    average_diff_fluxcom + base_theme,
    ncol = 2
  ) &
    theme(
      legend.position = 'bottom',
      legend.margin = margin(t = -0.4, r = 0, b = 0, l = 0, unit = 'cm')
    ),
  nrow = 2
) &
  theme(
    legend.title = element_text(
      size = 8,
      colour = '#23373b',
      margin = margin(0, 0, 0.05, 0, unit = 'cm')
    )
  )

ggsave_base(
  args$output,
  output,
  width = DISPLAY_SETTINGS$supplement_full_width,
  height = 11.5
)
