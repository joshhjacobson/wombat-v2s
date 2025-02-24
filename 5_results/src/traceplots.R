library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(patchwork)
library(ggrepel)

source(Sys.getenv('UTILS_PARTIAL'))
source(Sys.getenv('DISPLAY_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--samples')
parser$add_argument('--output')
args <- parser$parse_args()

samples <- readRDS(args$samples)

plot_traces <- function(x, names, title_shift) {
  if (!is.matrix(x)) x <- t(t(x))
  iterations <- (N_MCMC_WARM_UP + 1) : N_MCMC_SAMPLES
  df <- data.frame(
    iteration = rep(iterations, ncol(x)),
    value = as.vector(x[iterations, ])
  )

  if (!missing(names)) {
    df$name <- rep(names, each = N_MCMC_SAMPLES - N_MCMC_WARM_UP)

    output <- ggplot(mapping = aes(iteration, value, colour = name))
  } else {
    output <- ggplot(mapping = aes(iteration, value))
  }

  output <- output +
    geom_line(data = df)

  if (!missing(names)) {
    label_iteration <- round((N_MCMC_WARM_UP + 1 + N_MCMC_SAMPLES) / 2)
    df <- data.frame(
      name = names,
      iteration = label_iteration,
      value = colMeans(x[iterations, ])
    )
    output <- output + if (length(names) > 2) {
      geom_label_repel(
        data = df,
        size = 3,
        mapping = aes(label = name),
        direction = 'x',
        seed = 1
      )
    } else {
      geom_label(
        data = df,
        size = 3,
        mapping = aes(label = name)
      )
    }
  }

  output <- output +
    guides(colour = 'none') +
    scale_y_continuous(expand = expansion(mult = 0.15)) +
    labs(x = 'Iteration', y = 'Value', colour = NULL) +
    theme(
        plot.margin = margin(t = 1, r = 3, b = 0, l = 0, unit = 'mm'),
        plot.title = element_text(margin = margin(b = -0.5, unit = 'mm'))
      )

  if (!missing(title_shift)) {
    output <- output + theme(
      plot.title = element_text(margin = margin(t = title_shift, unit = 'mm'))
    )
  }

  output
}

to_bio_labels <- function(x) {
  c('bio_assim' = 'gpp', 'bio_resp_tot' = 'resp')[colnames(x)]
}

scale_bio_inventory <- scale_colour_manual(
  values = c('#018571', '#a6611a')
)

non_sif_groups <- colnames(samples$gamma)[!grepl('SIF', colnames(samples$gamma))]

output <- wrap_plots(
  plot_traces(
    samples$gamma[, non_sif_groups],
    c(
      '1_LNLG' = 'OCO-2 XCO2',
      'shipboard' = 'Shipboard',
      'tower' = 'Tower',
      'surface' = 'Surface',
      'aircraft' = 'Aircraft'
    )[non_sif_groups]
  ) +
    scale_colour_brewer(palette = 'Dark2') +
    ggtitle(expression(gamma[g]^Z)),
  plot_traces(samples$gamma[, '3_SIF']) +
    ggtitle(expression(gamma[sif]^Z)),
  plot_traces(samples$w_bio_season, to_bio_labels(samples$w_bio_season)) +
    ggtitle(expression(tau[c]^beta)) +
    scale_bio_inventory,
  plot_traces(samples$rho_bio_season) +
    ggtitle(expression(rho['gpp,resp']^beta)),
  plot_traces(samples$w_bio_resid, to_bio_labels(samples$w_bio_resid)) +
    ggtitle(expression(tau[c]^epsilon)) +
    scale_bio_inventory,
  plot_traces(samples$rho_bio_resid) +
    ggtitle(expression(rho['gpp,resp']^epsilon)),
  plot_traces(
    samples$alpha[, c('bio_assim.residual.Region05.2015-01', 'bio_resp_tot.residual.Region05.2015-01')],
    c('gpp', 'resp'),
    title_shift = 2
  ) +
    ggtitle(expression(alpha['c,6,5,5'])) +
    scale_bio_inventory,
  plot_traces(samples$kappa_bio_resid) +
    ggtitle(expression(kappa['bio']^epsilon)),
  plot_traces(
    samples$alpha[, c('bio_assim.intercept.Region02', 'bio_resp_tot.intercept.Region02')],
    c('gpp', 'resp')
  ) +
    ggtitle(expression(alpha['c,0,2'])) +
    scale_bio_inventory,
  plot_traces(
    samples$alpha[, c('bio_assim.trend.Region06', 'bio_resp_tot.trend.Region06')],
    c('gpp', 'resp')
  ) +
    ggtitle(expression(alpha['c,1,6'])) +
    scale_bio_inventory,
  ncol = 2
)

ggsave_base(
  args$output,
  output,
  width = DISPLAY_SETTINGS$supplement_full_width,
  height = 19
)
