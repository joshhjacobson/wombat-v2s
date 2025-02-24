library(ggplot2)

theme_set(theme_bw())
theme_replace(
  legend.background = element_blank(),
  legend.key = element_blank(),
  panel.background = element_blank(),
  strip.background = element_blank(),
  plot.background = element_blank(),
  # panel.border = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.minor.y = element_blank(),
  axis.text = element_text(colour = '#23373b'),
  axis.title = element_text(colour = '#23373b'),
  plot.margin = margin(t = 1, r = 0, b = 0, l = 1, unit = 'mm')
)

DISPLAY_SETTINGS <- list(
  full_width = 14.32,
  supplement_full_width = 16.5,
  full_height = 20,
  png_plot_dpi = 320,
  colour_key = c(
    'Truth' = 'black',
    'FLUXCOM' = '#e69f00cc',
    'Bottom-up' = '#8c8c8c',
    'Prior' = '#aa4499',
    'v2.0 posterior' = '#56b4e9cc',
    'v2.S posterior' = '#009e73',
    'Without SIF, fixed RLT' = '#d55e00',
    'With SIF, fixed RLT' = '#009e73',
    'Without SIF, inferred RLT' = '#d55e00',
    'With SIF, inferred RLT' = '#009e73',
    'Without SIF' = '#d55e00',
    'With SIF' = '#009e73'
  ),
  linetype_key = c(
    'Truth' = 'solid',
    'FLUXCOM' = '1131',
    'Bottom-up' = '11',
    'Prior' = '41',
    'v2.0 posterior' = '41',
    'v2.S posterior' = '41',
    'Without SIF, fixed RLT' = '41',
    'With SIF, fixed RLT' = '41',
    'Without SIF, inferred RLT' = '41',
    'With SIF, inferred RLT' = '41',
    'Without SIF' = '41',
    'With SIF' = '41'
  )
)

REGION_PLOT_SETTINGS <- list(
  'global' = list(
    latitude_lower = -90,
    latitude_upper = 90,
    lowercase_title = 'global',
    titlecase_title = 'Global',
    numeric_title = 'Global',
    in_supplement = FALSE
  ),
  'n-boreal' = list(
    latitude_lower = 50,
    latitude_upper = 90,
    lowercase_title = 'northern boreal (50°N-90°N)',
    titlecase_title = 'Northern boreal (50°N-90°N)',
    numeric_title = '50°N-90°N',
    in_supplement = TRUE
  ),
  'n-temperate' = list(
    latitude_lower = 23,
    latitude_upper = 50,
    lowercase_title = 'northern temperate (23°N-50°N)',
    titlecase_title = 'Northern temperate (23°N-50°N)',
    numeric_title = '23°N-50°N',
    in_supplement = TRUE
  ),
  'tropical' = list(
    latitude_lower = -23,
    latitude_upper = 23,
    lowercase_title = 'tropical (23°S-23°N)',
    titlecase_title = 'Tropical (23°S-23°N)',
    numeric_title = '23°S-23°N',
    in_supplement = TRUE
  ),
  'n-tropical' = list(
    latitude_lower = 0,
    latitude_upper = 23,
    lowercase_title = 'northern tropical (0°-23°N)',
    titlecase_title = 'Northern tropical (0°-23°N)',
    numeric_title = '0°-23°N',
    in_supplement = TRUE
  ),
  's-tropical' = list(
    latitude_lower = -23,
    latitude_upper = 0,
    lowercase_title = 'southern tropical (23°S-0°)',
    titlecase_title = 'Southern tropical (23°S-0°)',
    numeric_title = '23°S-0°',
    in_supplement = TRUE
  ),
  's-extratropical' = list(
    latitude_lower = -90,
    latitude_upper = -23,
    lowercase_title = 'southern extratropical (90°S-23°S)',
    titlecase_title = 'Southern extratropical (90°S-23°S)',
    numeric_title = '90°S-23°S',
    in_supplement = TRUE
  )
)

get_legend <- function(plot_object){
  tmp <- ggplot_gtable(ggplot_build(plot_object))
  legend_index <- which(sapply(tmp$grobs, function(x) x$name) == 'guide-box')
  tmp$grobs[[legend_index]]
}

ggsave_base <- function(filename, plot, bg = 'transparent', dpi = DISPLAY_SETTINGS$png_plot_dpi, ...) {
  ggsave(
    filename,
    plot,
    units = 'cm',
    dpi = dpi,
    bg = bg,
    ...
  )
}

ggsave_fullwidth <- function(...) {
  ggsave_base(..., width = DISPLAY_SETTINGS$full_width)
}

scale_fill_binned_custom <- function(
  palette_name,
  symmetric = FALSE,
  reverse = FALSE,
  trim_colours = FALSE,
  na.value = 'black',
  guide = 'coloursteps',
  aesthetics = 'fill',
  ...
) {
  n_colours <- ifelse(symmetric, 10, 9)
  palette_colours <- scico::scico(
    n_colours,
    palette = palette_name,
    direction = ifelse(reverse, -1, 1)
  )
  if (trim_colours) {
    palette_colours <- palette_colours[4:n_colours - 2]
  }
  palette <- scales::gradient_n_pal(palette_colours, NULL, 'Lab')
  binned_scale(
    aesthetics = 'fill',
    palette = function(x) {
      palette(seq(0, 1, length.out = length(x)))
    },
    na.value = na.value,
    guide = guide,
    ...
  )
}

plot_map <- function(
  df,
  variable,
  breaks,
  limits,
  palette,
  show_excess = TRUE,
  label_precision = 0,
  drop_second_labels = FALSE,
  symmetric = FALSE,
  reverse = FALSE,
  trim_colours = FALSE,
  bar_width = 11
) {
  base_labels <- sprintf(paste0('%.', label_precision, 'f'), breaks)
  labels <- if (show_excess) {
    ifelse(
      abs(breaks) == max(abs(breaks)),
      sprintf(
        paste0('%s%.', label_precision, 'f'),
        ifelse(breaks < 0, '< -', '> '),
        max(abs(breaks))
      ),
      base_labels
    )
  } else {
    base_labels
  }
  if (drop_second_labels) {
    labels[seq(2, length(breaks), by = 2)] <- ''
  }
  ggplot() +
    geom_sf(data = earth_bbox_sf, fill = '#dddddd', colour = 'black') +
    geom_stars(
      data = df,
      aes(fill = {{ variable }})
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
        x = c(-173, -180, -180),
        y = c(-23, 23, 50),
        label = c('23°S', '23°N', '50°N')
      ),
      mapping = aes(x = x, y = y, label = label),
      size = 7/.pt,
      nudge_y = 6
    ) +
    coord_sf(
      crs = sf::st_crs('ESRI:54012'),
      default_crs = sf::st_crs('WGS84')
    ) +
    scale_fill_binned_custom(
      palette,
      breaks = breaks,
      limits = limits,
      labels = labels,
      symmetric = symmetric,
      reverse = reverse,
      trim_colours = trim_colours,
      guide = guide_coloursteps(
        title.position = 'top',
        title.hjust = 0.5,
        label.theme = element_text(
          size = 7,
          colour = '#23373b',
          margin = margin(0.1, 0, 0, 0, unit = 'cm')
        ),
        frame.colour = NA,
        barwidth = bar_width,
        barheight = 0.5,
        even.steps = TRUE
      ),
      na.value = NA
    ) +
    theme(
      panel.border = element_blank()
    ) +
    labs(x = NULL, y = NULL)
}
