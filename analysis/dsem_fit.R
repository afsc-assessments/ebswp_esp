# dsem_fit.R
# Confirmatory DSEM estimation for the EBS pollock early-survival DAGs.
# Companion to analysis/pcmci_plus.py: PCMCI+ screens links, DSEM estimates
# effect sizes on each specified DAG and ranks the four hypotheses by AIC.
# See dsem.qmd for the narrative. Reference: Thorson et al. (2024); Ma et al. (2026).

library(dsem)
library(dplyr)
library(here)

# -----------------------------------------------------------------------------
# 1. Cohort table (shared input with the PCMCI+ workflow)
# -----------------------------------------------------------------------------
# Lags are encoded in the column names (cohort alignment):
#   *_t  = age-0 / spawning year, *_t1 = age-1 year, *_t3 = age-3 recruitment.
cohort <- readr::read_csv(here("data", "cohort_table.csv"), show_col_types = FALSE)

# Standardize to mean 0 / sd 1 (comparable slopes); retain NAs for the
# state-space treatment in dsem.
z     <- function(x) as.numeric(scale(x))
nodes <- setdiff(names(cohort), "year")
tsdat <- ts(as.data.frame(lapply(cohort[nodes], z)), start = min(cohort$year))

# -----------------------------------------------------------------------------
# 2. The four DAGs as DSEM path specifications ('from -> to, lag, name')
# -----------------------------------------------------------------------------
dag_A <- "
  SSB                  -> age3_recruitment_t3,  0, ssb_R
  ice_retreat_t        -> copepod_euphausiid_t, 0, ice_cope
  wind_mixing_t        -> copepod_euphausiid_t, 0, wind_cope
  late_summer_SST_t    -> copepod_euphausiid_t, 0, sst_cope
  copepod_euphausiid_t -> age0_energy_t,        0, cope_en
  age0_energy_t        -> age3_recruitment_t3,  0, en_R
"

dag_B <- "
  SSB                          -> age3_recruitment_t3,         0, ssb_R
  adult_pollock_biomass_t1     -> juvenile_predator_overlap_t1, 0, adult_ov
  arrowtooth_biomass_t1        -> juvenile_predator_overlap_t1, 0, atf_ov
  cold_pool_t1                 -> juvenile_predator_overlap_t1, 0, cp_ov
  late_summer_SST_t            -> juvenile_predator_overlap_t1, 0, sst_ov
  juvenile_predator_overlap_t1 -> age3_recruitment_t3,          0, ov_R
"

# DAG C requires the transport columns (spawning_location_t, larval_transport_t,
# prey_match_t, age0_abundance_t). Enable once they exist in the cohort table.
dag_C <- "
  SSB                 -> age3_recruitment_t3, 0, ssb_R
  late_summer_SST_t   -> spawning_location_t, 0, sst_spawn
  wind_mixing_t       -> larval_transport_t,  0, wind_tran
  spawning_location_t -> larval_transport_t,  0, spawn_tran
  larval_transport_t  -> prey_match_t,        0, tran_prey
  larval_transport_t  -> age0_abundance_t,    0, tran_a0
  prey_match_t        -> age0_abundance_t,    0, prey_a0
  age0_abundance_t    -> age3_recruitment_t3, 0, a0_R
"

dag_D <- "
  ice_retreat_t                -> copepod_euphausiid_t,         0, ice_cope
  copepod_euphausiid_t         -> age0_energy_t,                0, cope_en
  age0_energy_t                -> age1_index_t1,                0, en_a1
  age1_index_t1                -> juvenile_predator_overlap_t1, 0, a1_ov
  adult_pollock_biomass_t1     -> juvenile_predator_overlap_t1, 0, adult_ov
  juvenile_predator_overlap_t1 -> age3_recruitment_t3,          0, ov_R
  age3_recruitment_t3          -> age3_recruitment_t3,          1, ar_R
"
# DAG D's switching claim is state-dependence; add an interaction column to the
# cohort table (e.g. ice_x_adult = z(ice_retreat_t * adult_pollock_biomass_t1))
# and the path 'ice_x_adult -> age3_recruitment_t3, 0, switch', or fit A vs B by
# warm/cold regime. The nonstationarity check in section 4 tests this continuously.

# -----------------------------------------------------------------------------
# 3. Fit and rank by AIC (the confirmatory step)
# -----------------------------------------------------------------------------
fit_dag <- function(sem, data = tsdat) {
  dsem(sem = sem, tsdata = data, estimate_delta0 = TRUE,
       family = rep("normal", ncol(data)))
}

fits <- list(
  A_bottom_up = fit_dag(dag_A),
  B_top_down  = fit_dag(dag_B),
  # C_transport = fit_dag(dag_C),   # enable when transport columns exist
  D_switching = fit_dag(dag_D)
)

aic_table <- tibble(
  dag  = names(fits),
  aic  = sapply(fits, AIC)
) |>
  mutate(d_aic = aic - min(aic)) |>
  arrange(aic)

best <- fits[[aic_table$dag[1]]]

# Persist outputs for the report.
dir.create(here("outputs"), showWarnings = FALSE)
readr::write_csv(aic_table, here("outputs", "dsem_aic.csv"))
# summary(best) returns the standardized path coefficients with SE and p:
readr::write_csv(as.data.frame(summary(best)), here("outputs", "dsem_paths.csv"))

print(aic_table)

# -----------------------------------------------------------------------------
# 4. Nonstationarity check (the highest-value sensitivity, after Ma et al.)
# -----------------------------------------------------------------------------
# Rolling-window slopes for a focal DAG; a sign change supports switching control.
rolling_slopes <- function(sem, win = 21, data = tsdat) {
  yrs <- as.numeric(time(data))
  out <- lapply(seq_len(length(yrs) - win + 1), function(i) {
    sub <- window(data, start = yrs[i], end = yrs[i + win - 1])
    ft  <- tryCatch(fit_dag(sem, sub), error = function(e) NULL)
    if (is.null(ft)) return(NULL)
    s <- as.data.frame(summary(ft))
    s$center_year <- yrs[i] + (win - 1) / 2
    s
  })
  dplyr::bind_rows(out)
}
# ns_B <- rolling_slopes(dag_B)
# readr::write_csv(ns_B, here("outputs", "dsem_nonstationarity_B.csv"))
