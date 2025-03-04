library(argparse)
library(dplyr, warn.conflicts = FALSE)

parser <- ArgumentParser()
parser$add_argument('--samples-wombat-v2')
parser$add_argument('--output')
args <- parser$parse_args()

samples_wombat_v2 <- readRDS(args$samples_wombat_v2)

output <- samples_wombat_v2$alpha_df %>%
  select(-value_samples)

fst::write_fst(output, args$output)
