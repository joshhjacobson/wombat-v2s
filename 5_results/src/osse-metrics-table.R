library(argparse)
library(dplyr, warn.conflicts = FALSE)

source(Sys.getenv('UTILS_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--metric')
parser$add_argument('--digits', type = 'integer', default = 2)
parser$add_argument('--flux-samples-alpha0-fixresp-wsif')
parser$add_argument('--flux-samples-alpha0-fixresp-wosif')
parser$add_argument('--flux-samples-alpha0-freeresp-wsif')
parser$add_argument('--flux-samples-alpha0-freeresp-wosif')
parser$add_argument('--flux-samples-alphav2-fixresp-wsif')
parser$add_argument('--flux-samples-alphav2-fixresp-wosif')
parser$add_argument('--flux-samples-alphav2-freeresp-wsif')
parser$add_argument('--flux-samples-alphav2-freeresp-wosif')
parser$add_argument('--flux-samples-alphap-fixresp-wsif')
parser$add_argument('--flux-samples-alphap-fixresp-wosif')
parser$add_argument('--flux-samples-alphap-freeresp-wsif')
parser$add_argument('--flux-samples-alphap-freeresp-wosif')
parser$add_argument('--flux-samples-alphan-fixresp-wsif')
parser$add_argument('--flux-samples-alphan-fixresp-wosif')
parser$add_argument('--flux-samples-alphan-freeresp-wsif')
parser$add_argument('--flux-samples-alphan-freeresp-wosif')
parser$add_argument('--output')
args <- parser$parse_args()

printf <- function(...) cat(sprintf(...))
collapse0 <- function(x) paste0(x, collapse = '')
paste_columns <- function(x) paste0(x, collapse = ' & ')
paste_columns_fmt <- function(x, is_bold, digits = args$digits) {
  format_digits <- sprintf(paste0('%.', digits, 'f'), x)
  format_bold <- ifelse(is_bold, paste0('\\textbf{', format_digits, '}'), format_digits)
  paste0(format_bold, collapse = ' & ')
}
read_flux_samples <- function(filename, estimates = 'Posterior') {
  readRDS(filename) %>% filter(estimate %in% estimates)
}

log_debug('Loading flux samples')
flux_aggregates_samples <- bind_rows(
  read_flux_samples(
    args$flux_samples_alpha0_fixresp_wsif,
    c('Truth', 'Posterior')
  ) %>%
    mutate(
      truth = 'Bottom-up',
      resp_term = 'Fixed',
      estimate = if_else(
        estimate == 'Posterior',
        'With SIF',
        estimate
      )
    ),
  read_flux_samples(
    args$flux_samples_alpha0_fixresp_wosif,
  ) %>%
    mutate(
      truth = 'Bottom-up',
      resp_term = 'Fixed',
      estimate = 'Without SIF'
    ),
  read_flux_samples(
    args$flux_samples_alphav2_fixresp_wsif,
    c('Truth', 'Posterior')
  ) %>%
    mutate(
      truth = 'v2.0 mean',
      resp_term = 'Fixed',
      estimate = if_else(
        estimate == 'Posterior',
        'With SIF',
        estimate
      )
    ),
  read_flux_samples(
    args$flux_samples_alphav2_fixresp_wosif,
  ) %>%
    mutate(
      truth = 'v2.0 mean',
      resp_term = 'Fixed',
      estimate = 'Without SIF'
    ),
  read_flux_samples(
    args$flux_samples_alphap_fixresp_wsif,
    c('Truth', 'Posterior')
  ) %>%
    mutate(
      truth = 'Positive shift',
      resp_term = 'Fixed',
      estimate = if_else(
        estimate == 'Posterior',
        'With SIF',
        estimate
      )
    ),
  read_flux_samples(
    args$flux_samples_alphap_fixresp_wosif,
  ) %>%
    mutate(
      truth = 'Positive shift',
      resp_term = 'Fixed',
      estimate = 'Without SIF'
    ),
  read_flux_samples(
    args$flux_samples_alphan_fixresp_wsif,
    c('Truth', 'Posterior')
  ) %>%
    mutate(
      truth = 'Negative shift',
      resp_term = 'Fixed',
      estimate = if_else(
        estimate == 'Posterior',
        'With SIF',
        estimate
      )
    ),
  read_flux_samples(
    args$flux_samples_alphan_fixresp_wosif,
  ) %>%
    mutate(
      truth = 'Negative shift',
      resp_term = 'Fixed',
      estimate = 'Without SIF'
    ),
  read_flux_samples(
    args$flux_samples_alpha0_freeresp_wsif,
    c('Truth', 'Posterior')
  ) %>%
    mutate(
      truth = 'Bottom-up',
      resp_term = 'Inferred',
      estimate = if_else(
        estimate == 'Posterior',
        'With SIF',
        estimate
      )
    ),
  read_flux_samples(
    args$flux_samples_alpha0_freeresp_wosif,
  ) %>%
    mutate(
      truth = 'Bottom-up',
      resp_term = 'Inferred',
      estimate = 'Without SIF'
    ),
  read_flux_samples(
    args$flux_samples_alphav2_freeresp_wsif,
    c('Truth', 'Posterior')
  ) %>%
    mutate(
      truth = 'v2.0 mean',
      resp_term = 'Inferred',
      estimate = if_else(
        estimate == 'Posterior',
        'With SIF',
        estimate
      )
    ),
  read_flux_samples(
    args$flux_samples_alphav2_freeresp_wosif,
  ) %>%
    mutate(
      truth = 'v2.0 mean',
      resp_term = 'Inferred',
      estimate = 'Without SIF'
    ),
  read_flux_samples(
    args$flux_samples_alphap_freeresp_wsif,
    c('Truth', 'Posterior')
  ) %>%
    mutate(
      truth = 'Positive shift',
      resp_term = 'Inferred',
      estimate = if_else(
        estimate == 'Posterior',
        'With SIF',
        estimate
      )
    ),
  read_flux_samples(
    args$flux_samples_alphap_freeresp_wosif,
  ) %>%
    mutate(
      truth = 'Positive shift',
      resp_term = 'Inferred',
      estimate = 'Without SIF'
    ),
  read_flux_samples(
    args$flux_samples_alphan_freeresp_wsif,
    c('Truth', 'Posterior')
  ) %>%
    mutate(
      truth = 'Negative shift',
      resp_term = 'Inferred',
      estimate = if_else(
        estimate == 'Posterior',
        'With SIF',
        estimate
      )
    ),
  read_flux_samples(
    args$flux_samples_alphan_freeresp_wosif,
  ) %>%
    mutate(
      truth = 'Negative shift',
      resp_term = 'Inferred',
      estimate = 'Without SIF'
    )
) %>%
  filter(inventory != 'Total')

log_debug('Computing metrics')
osse_fluxes <- flux_aggregates_samples %>%
  filter(estimate != 'Truth') %>%
  left_join(
    flux_aggregates_samples %>%
      filter(estimate == 'Truth') %>%
      distinct(inventory, truth, region, time, flux_truth = flux_mean),
    by = c('inventory', 'truth', 'region', 'time')
  ) %>%
  mutate(
    truth = factor(
      truth,
      levels = c('Bottom-up', 'v2.0 mean', 'Positive shift', 'Negative shift')
    ),
    rlt_switch = factor(
      resp_term,
      levels = c('Fixed', 'Inferred')
    ),
    sif_switch = factor(
      estimate,
      levels = c('With SIF', 'Without SIF')
    ),
    inventory = factor(c(
      'GPP' = 'GPP',
      'Respiration' = 'Resp.',
      'NEE' = 'NEE',
      'Ocean' = 'Ocean'
    )[inventory], levels = c('GPP', 'Resp.', 'NEE', 'Ocean'))
  )

compute_metric <- function(metric, mean, samples, truth) {
  switch(
    metric,
    rmse = sqrt(mean((mean - truth)^2)),
    crps = mean(scoringRules::crps_sample(truth, samples))
  )
}

metrics <- osse_fluxes %>%
  group_by(truth, rlt_switch, sif_switch, inventory) %>%
  summarise(
    value = compute_metric(
      args$metric,
      flux_mean,
      flux_samples,
      flux_truth
    ),
    .groups = 'drop'
  ) %>%
  arrange(truth, rlt_switch, sif_switch, inventory) %>%
  tidyr::pivot_wider(
    names_from = inventory,
    values_from = value
  )

truth_cases <- levels(metrics$truth)
rlt_cases <- levels(metrics$rlt_switch)
sif_cases <- levels(metrics$sif_switch)
inventories <-  levels(osse_fluxes$inventory)

is_group_min <- function(x) {
  x_rounded <- round(x, args$digits)
  output <- rep(FALSE, length(x))
  if (sum(x_rounded == min(x_rounded)) == 1) {
    output[x_rounded == min(x_rounded)] <- TRUE
  }
  output
}

format_bold <- metrics %>%
  group_by(truth) %>%
  mutate(across(all_of(inventories), is_group_min)) %>%
  ungroup()

format_truth_case <- function(case) {
  if (case %in% c('Bottom-up', 'v2.0 mean')) {
    paste0(as.character(case), '\\textsuperscript{b}')
  } else {
    as.character(case)
  }
}

metric_name <- function(metric) {
  stopifnot(tolower(metric) %in% c('rmse', 'crps'))
  switch(
    tolower(metric),
    rmse = 'RMSE',
    crps = 'Mean CRPS'
  )
}

log_debug('Writing table to {args$output}')
sink(args$output)
printf('\\begin{tabular}{lll%s}\n\\toprule\n', collapse0(rep('C', length(inventories))))
printf(
  'True Flux & \\multicolumn{2}{c}{Inversion Setup} & \\multicolumn{%d}{c}{%s [PgC/year]} \\\\\n',
  length(inventories),
  metric_name(args$metric)
)
printf('\\cmidrule(lr){2-3} \\cmidrule(l{10pt}){4-%d}\n', ncol(metrics))
printf('& RLT\\textsuperscript{a} & Includes SIF & %s \\\\\n\\midrule\n', paste_columns(inventories))
for (i in seq_len(nrow(metrics))) {
  metrics_row <- metrics[i, ]
  format_row <- format_bold[i, ]
  printf(
    '%s & %s & %s & %s \\\\\n',
    ifelse(
      metrics_row$rlt_switch == rlt_cases[1] & metrics_row$sif_switch == sif_cases[1],
      format_truth_case(metrics_row$truth),
      ''
    ),
    ifelse(metrics_row$sif_switch == sif_cases[1], as.character(metrics_row$rlt_switch), ''),
    ifelse(metrics_row$sif_switch == 'With SIF', 'Yes', 'No'),
    paste_columns_fmt(metrics_row[inventories], format_row[inventories])
  )
  if (i != nrow(metrics)) {
    if (metrics_row$truth != metrics[i + 1, ]$truth) cat('\\midrule\n')
  }
}
cat('\\bottomrule\n\\end{tabular}\n')
sink(NULL)

log_debug('Done')
