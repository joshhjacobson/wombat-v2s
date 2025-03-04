library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(Rcpp)

rcpp_cache_dir <- Sys.getenv('RCPP_CACHE_DIR')
options(rcpp.cache.dir = if (rcpp_cache_dir == '') tempdir() else rcpp_cache_dir)

source(Sys.getenv('UTILS_PARTIAL'))
sourceCpp(Sys.getenv('UTILS_CPP_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--region-mask')
parser$add_argument('--sib4-climatology-assim')
parser$add_argument('--sib4-climatology-resp-tot')
parser$add_argument('--basis-vectors')
parser$add_argument('--control-emissions')
parser$add_argument('--alpha-wombat-v2')
parser$add_argument('--fix-resp-linear')
parser$add_argument('--delta', type = 'numeric')
parser$add_argument('--output')
args <- parser$parse_args()

basis_vectors <- fst::read_fst(args$basis_vectors)
control_emissions <- fst::read_fst(args$control_emissions)
alpha <- fst::read_fst(args$alpha_wombat_v2)

regions_land <- sprintf('Region%02d', 1:11)
regions_free <- setdiff(regions_land, args$fix_resp_linear)
stopifnot(length(regions_free) > 0)

region_list <- readRDS(args$region_mask)
region_values <- region_list[[1]]
for (i in 2:length(region_list)) {
  x <- region_list[[i]]
  region_values[x == 1] <- i
}
region_values[region_values == 0] <- NA

region_grid <- expand.grid(
  longitude = attr(region_list$Region01, 'longitude'),
  latitude = attr(region_list$Region01, 'latitude')
) %>%
  mutate(
    region = as.vector(region_values)
  ) %>%
  filter(!is.na(region)) %>%
  mutate(
    region = factor(names(region_list)[region])
  )

climatology_linear <- bind_rows(
  read_climatology(args$sib4_climatology_assim) %>%
    mutate(
      value = -value,
      inventory = factor('bio_assim')
    ),
  read_climatology(args$sib4_climatology_resp_tot) %>%
    mutate(inventory = factor('bio_resp_tot'))
) %>%
  rename(component = variable) %>%
  filter(component %in% c('intercept', 'trend')) %>%
  left_join(
    region_grid,
    by = c('longitude', 'latitude')
  ) %>%
  filter(region %in% regions_free) %>%
  left_join(
    control_emissions %>% distinct(longitude, latitude, area),
    by = c('longitude', 'latitude')
  ) %>%
  mutate(month = factor(NA)) %>%
  add_basis_vector(basis_vectors)

perturbations <- climatology_linear %>%
  group_by(basis_vector) %>%
  summarise(value = sum(area * value)) %>%
  left_join(
    climatology_linear %>%
      distinct(basis_vector, inventory, component, region),
    by = 'basis_vector'
  )

log_debug('Computing adjusted linear component for bio_assim')
alpha_bio_assim_linear <- alpha %>%
  filter(
    inventory == 'bio_assim',
    component %in% c('intercept', 'trend'),
    region %in% regions_free
  ) %>%
  mutate(
    delta = args$delta,
    value = value + delta
  ) %>%
  select(c(basis_vector, inventory, component, region, value, delta))

log_debug('Computing adjusted linear component for bio_resp_tot')
# NOTE(jhj): This construction ensures that the implied linear component for NEE
# is unchanged at regional scales or larger
alpha_bio_resp_tot_linear <- perturbations %>%
  left_join(
    alpha_bio_assim_linear %>% select(
      c(basis_vector, delta)
    ),
    by = 'basis_vector'
  ) %>%
  mutate(
    value = if_else(inventory == 'bio_assim', delta * value, value)
  ) %>%
  select(-c(basis_vector, delta)) %>%
  tidyr::pivot_wider(
    names_from = inventory,
    values_from = value
  ) %>%
  mutate(
    value = -bio_assim / bio_resp_tot,
    inventory = 'bio_resp_tot'
  ) %>%
  select(-c(bio_assim, bio_resp_tot)) %>%
  left_join(
    perturbations %>%
      distinct(basis_vector, inventory, component, region),
    by = c('inventory', 'component', 'region')
  )

alpha_bio_linear <- bind_rows(
  alpha_bio_assim_linear %>% select(-delta),
  alpha_bio_resp_tot_linear
) %>%
  mutate(
    basis_vector_str = as.character(basis_vector),
    month = NA
  )

output <- bind_rows(
  alpha_bio_linear,
  alpha %>%
    filter(
      !(
        inventory == 'bio_assim' &
          component %in% c('intercept', 'trend') &
          region %in% regions_free
      )
    )
) %>%
  mutate(
    basis_vector = factor(
      basis_vector_str,
      levels = levels(basis_vectors$basis_vector)
    )
  )

log_debug('Writing to {args$output}')
fst::write_fst(output, args$output)

log_debug('Done')
