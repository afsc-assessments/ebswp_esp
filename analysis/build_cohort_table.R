# build_cohort_table.R
# Assemble data/cohort_table.csv: one row per year class (cohort), every covariate
# aligned to the cohort's birth (age-0) year. Shared input for the PCMCI+
# (pcmci.qmd) and DSEM (dsem.qmd) workflows. See data/cohort_table_dictionary.qmd.
#
# Status: the ebswp + coldpool columns below are RUNNABLE against the accepted run.
# The ESP/BASIS columns (age-0 energy, juvenile-predator overlap, ice/wind,
# arrowtooth, age-1 index) are left NA until those exports are located.
#
# Conventions (AGENTS.md): 2-space indent, snake_case, here() for paths.

library(dplyr)
library(here)
library(ebswp)      # provides read_admb()
# library(coldpool) # loaded in section 2

# Accepted assessment run (ADMB 'pm' model). read_admb() wants the path with no
# extension; it reads pm.rep + pm.std.
accepted_run <- "/Users/jim/_mymods/pollock/ebswp_assessment/2025/runs/base/pm"

# Helper: pull a year-indexed series from an sdreport-style matrix
# (col 1 = year, col 2 = estimate) into a tidy 2-column frame.
col2 <- function(mat, value_name) {
  tibble(year = mat[, 1], !!value_name := mat[, 2])
}

# -----------------------------------------------------------------------------
# 1. ebswp assessment series (SSB, recruitment, adult abundance, copepod index)
# -----------------------------------------------------------------------------
A <- read_admb(accepted_run)

ssb   <- col2(A$SSB,    "ssb")        # spawning stock biomass, by calendar year
nage3 <- col2(A$Nage_3, "nage3")      # numbers at age 3, by calendar year
# Copepod index is carried INSIDE the assessment (fit covariate), years 2002+:
cope  <- tibble(year = A$yrs_cope, copepod_euphausiid = A$obs_cope)
# age-1 recruitment is also available if you prefer it as the response:
rec1  <- col2(A$R, "age1_recruit")    # optional; see note in section 4

assess <- ssb |>
  full_join(nage3, by = "year") |>
  full_join(cope,  by = "year") |>
  arrange(year)

# -----------------------------------------------------------------------------
# 2. coldpool package (Rohan et al. 2022): cold pool extent + temperatures
# -----------------------------------------------------------------------------
# install.packages("coldpool",
#   repos = c("https://afsc-gap-products.r-universe.dev", getOption("repos")))
# Verify column names with: names(coldpool::cold_pool_index)
# Prefer the installed package; otherwise use the bundled snapshot in data-raw/
# (cold_pool_index.csv, pulled from github.com/afsc-gap-products/coldpool).
cpi_local <- here("data-raw", "cold_pool_index.csv")
if (requireNamespace("coldpool", quietly = TRUE)) {
  cpi <- as.data.frame(coldpool::cold_pool_index)
} else if (file.exists(cpi_local)) {
  message("Using bundled data-raw/cold_pool_index.csv (coldpool not installed).")
  cpi <- readr::read_csv(cpi_local, show_col_types = FALSE)
} else {
  stop("No cold pool source: install coldpool or provide data-raw/cold_pool_index.csv")
}
phys <- cpi |>
  transmute(year            = YEAR,
            cold_pool        = AREA_LTE2_KM2,            # bottom area <= 2 C
            late_summer_SST  = MEAN_SURFACE_TEMPERATURE,
            bottom_temp      = MEAN_GEAR_TEMPERATURE)

# -----------------------------------------------------------------------------
# 3. ESP / BASIS gaps (not on this machine yet) -- tidy CSVs with a 'year' col
# -----------------------------------------------------------------------------
# Each should be year-indexed; join then shift in section 4. Until they exist,
# they enter as NA so DAG A can run on the ebswp/coldpool columns alone.
#   age0_energy           Heintz et al. 2013; Moss et al. 2009 (BASIS fall energy)
#   ice_retreat, wind_mixing   ESR physical indicators
#   age1_index            age-1 survey index
#   arrowtooth_biomass    arrowtooth assessment/survey
#   juvenile_predator_overlap  VAST overlap product (project's key build)
gap_path <- here("data-raw")        # drop exports here when available
read_gap <- function(file, value) {
  f <- file.path(gap_path, file)
  if (file.exists(f)) readr::read_csv(f, show_col_types = FALSE)
  else tibble(year = integer(), !!value := numeric())
}
age0_energy <- read_gap("age0_energy.csv",        "age0_energy")
ice         <- read_gap("ice_retreat.csv",        "ice_retreat")
wind        <- read_gap("wind_mixing.csv",        "wind_mixing")
age1_idx    <- read_gap("age1_index.csv",         "age1_index")
arrowtooth  <- read_gap("arrowtooth_biomass.csv", "arrowtooth_biomass")
overlap     <- read_gap("predator_overlap.csv",   "juvenile_predator_overlap")

# -----------------------------------------------------------------------------
# 4. Join on calendar year, then shift to cohort alignment and write
# -----------------------------------------------------------------------------
cal <- assess |>
  full_join(phys,        by = "year") |>
  left_join(age0_energy, by = "year") |>
  left_join(ice,         by = "year") |>
  left_join(wind,        by = "year") |>
  left_join(age1_idx,    by = "year") |>
  left_join(arrowtooth,  by = "year") |>
  left_join(overlap,     by = "year") |>
  arrange(year)

# Cohort alignment. Row `year` = cohort's age-0 year.
#   *_t  : same calendar year   -> take as-is
#   *_t1 : cohort age-1 year    -> lead(x, 1)
#   *_t3 : cohort age-3 year    -> lead(x, 3)
# Adult predator node uses SSB at age-1 year (mature biomass = cannibalism proxy);
# response uses numbers-at-age-3 three years on (the cohort's own age-3).
cohort <- cal |>
  transmute(
    year,
    SSB                          = ssb,
    ice_retreat_t                = ice_retreat,
    late_summer_SST_t            = late_summer_SST,
    wind_mixing_t                = wind_mixing,
    copepod_euphausiid_t         = copepod_euphausiid,
    age0_energy_t                = age0_energy,
    cold_pool_t1                 = lead(cold_pool, 1),
    age1_index_t1                = lead(age1_index, 1),
    adult_pollock_biomass_t1     = lead(ssb, 1),          # alt: lead(nage3, 1)
    arrowtooth_biomass_t1        = lead(arrowtooth_biomass, 1),
    juvenile_predator_overlap_t1 = lead(juvenile_predator_overlap, 1),
    age3_recruitment_t3          = lead(nage3, 3)         # alt response: age-1 R
  ) |>
  filter(year >= 1977)

dir.create(here("data"), showWarnings = FALSE)
readr::write_csv(cohort, here("data", "cohort_table.csv"))
message("Wrote ", nrow(cohort), " cohorts to data/cohort_table.csv")
