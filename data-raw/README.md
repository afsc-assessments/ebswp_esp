# Source data

This directory holds source snapshots used to build the cohort-aligned analysis
table. Keep source field names and units intact where practical; transformations
belong in `analysis/build_cohort_table.R`.

## Current source

| File | Source | Role |
|---|---|---|
| `cold_pool_index.csv` | Snapshot of `coldpool::cold_pool_index` from the [AFSC coldpool project](https://github.com/afsc-gap-products/coldpool) | Fallback source for cold-pool area and surface and bottom temperature |

The build script preferentially uses an installed `coldpool` package and falls
back to this snapshot. Record the retrieval date or package version here when
the snapshot is refreshed.

## Expected indicator inputs

Optional files use one row per calendar year and the two columns `year` and the
value name shown below:

| File | Value column |
|---|---|
| `age0_energy.csv` | `age0_energy` |
| `ice_retreat.csv` | `ice_retreat` |
| `wind_mixing.csv` | `wind_mixing` |
| `age1_index.csv` | `age1_index` |
| `arrowtooth_biomass.csv` | `arrowtooth_biomass` |
| `predator_overlap.csv` | `juvenile_predator_overlap` |

For each added file, document its provider, retrieval date, units, spatial and
seasonal scope, missing-value convention, and any restrictions on redistribution.
