# Analysis-ready data

This directory contains cohort-aligned inputs used by PCMCI+ and DSEM.

- `cohort_table.csv` is the current assembled input. One row is one birth-year
  cohort; suffixes such as `_t1` and `_t3` identify the calendar year from which
  a value is read.
- `cohort_table_template.csv` is the column template for assembling or checking
  an input by hand.
- `cohort_table_dictionary.qmd` is the authoritative column contract, including
  timing, candidate sources, and optional transport variables.

Regenerate `cohort_table.csv` with:

```bash
Rscript analysis/build_cohort_table.R
```

Do not edit the generated table by hand. Update a source in `data-raw/` or the
transformation script, then rebuild so the change remains traceable.
