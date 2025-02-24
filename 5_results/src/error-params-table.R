library(argparse)
library(dplyr, warn.conflicts = FALSE)

parser <- ArgumentParser()
parser$add_argument('--hyperparameter-estimates')
parser$add_argument('--output')
args <- parser$parse_args()

ell_unit_key <- data.frame(
  hyperparameter_group = c(
    'Aircraft',
    'Shipboard',
    'Surface',
    'Tower',
    'OCO-2 XCO2',
    'OCO-2 SIF'
  ),
  scale_factor = c(
    24 * 60,
    24,
    24,
    24,
    60,
    60
  ),
  ell_unit = c(
    'minutes',
    'hours',
    'hours',
    'hours',
    'seconds',
    'seconds'
  )
)

hyperparameter_estimates <- fst::read_fst(args$hyperparameter_estimates) %>%
  filter(hyperparameter_group != 'lauder') %>%
  select(-ell_unit) %>%
  mutate(
    hyperparameter_group = factor(c(
      'aircraft' = 'Aircraft',
      'shipboard' = 'Shipboard',
      'surface' = 'Surface',
      'tower' = 'Tower',
      '1_LNLG' = 'OCO-2 XCO2',
      '3_SIF' = 'OCO-2 SIF'
    )[as.character(hyperparameter_group)], levels = c(
      'Aircraft',
      'Shipboard',
      'Surface',
      'Tower',
      'OCO-2 XCO2',
      'OCO-2 SIF'
    ))
  ) %>%
  arrange(hyperparameter_group) %>%
  mutate(hyperparameter_group = as.character(hyperparameter_group)) %>%
  left_join(ell_unit_key, by = 'hyperparameter_group') %>%
  mutate(ell = ell * scale_factor)


sink(args$output)
cat('\\begin{tabular}{lccl}\n\\toprule\n')
cat('Observation group & $\\hat{\\gamma}^Z_g$ & $\\hat{\\rho}^Z_g$ & \\multicolumn{1}{c}{$\\hat{\\ell}^Z_g$} \\\\\n\\midrule\n')
for (i in seq_len(nrow(hyperparameter_estimates))) {
  cat(sprintf(
    '%s & %.03f & %.03f & %.01f %s \\\\\n',
    sprintf(
      '%s',
      ifelse(
        hyperparameter_estimates$hyperparameter_group[i] == 'OCO-2 XCO2',
        'OCO-2 XCO\\textsubscript{2}',
        hyperparameter_estimates$hyperparameter_group[i]
      )
    ),
    hyperparameter_estimates$gamma[i],
    hyperparameter_estimates$rho[i],
    hyperparameter_estimates$ell[i],
    hyperparameter_estimates$ell_unit[i]
  ))
}
cat('\\bottomrule\n\\end{tabular}\n')
sink(NULL)
