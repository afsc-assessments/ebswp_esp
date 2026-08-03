# Analysis workflows

The analysis scripts share `data/cohort_table.csv`, whose timing convention is
defined in `data/cohort_table_dictionary.qmd`.

## Entry points

1. `build_cohort_table.R` joins assessment, cold-pool, and optional ESP/BASIS
   series by calendar year, then aligns them to birth-year cohorts.
2. `dsem_fit.R` standardizes the cohort series, fits the candidate confirmatory
   dynamic structural equation models, compares them with AIC, and writes model
   summaries to `outputs/`.

Run scripts from the repository root:

```bash
Rscript analysis/build_cohort_table.R
Rscript analysis/dsem_fit.R
```

The build currently depends on the local `ebswp` package and accepted assessment
run identified near the top of `build_cohort_table.R`. The DSEM fit requires the
`dsem`, `dplyr`, `here`, and `readr` packages. Several candidate indicator series
remain explicit gaps; see `data-raw/README.md` and the website's DSEM page.

Generated model outputs belong in `outputs/` and should not be hand-edited.
