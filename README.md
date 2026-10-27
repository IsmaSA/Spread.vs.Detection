# Separating spread from detection in freshwater invasion dynamics across Europe

R code to reproduce the analyses, tables and figures of the paper.

## Requirements

R 4.5.0 with data.table, sf, units, rgbif, CoordinateCleaner, rnaturalearth, rnaturalearthdata, sandwich, ggplot2, patchwork, scales, ragg and svglite. Versions used are listed in `session_info.txt`.

## Data

Occurrences come from 24 GBIF downloads (14 focal species and 10 target groups) listed with DOI and download key in `data/gbif_downloads.csv`. `R/01_data.R` retrieves them (about 1.7 GB), cleans the records and writes one file per taxon to `data/occ/`. `data/species.csv` gives the target group, year of introduction in Europe and whether native-range records were removed (native_box = 1) for each species.

## Running

From this folder:

```
Rscript run_all.R
```

Steps 02 to 04 can be split over several processes, e.g. `Rscript R/02_rates.R 1 4`, ..., `Rscript R/02_rates.R 4 4`.

| Script | Output |
|---|---|
| `R/01_data.R` | `data/occ/`: cleaned records in EPSG:3035 |
| `R/02_rates.R` | `output/runs/`: range curves and spread rates, 14 species x 10, 25 and 50 km grids |
| `R/03_radius.R` | `output/runs/`: accessible-area rarefaction with 75 and 300 km radii |
| `R/04_simulation.R` | `output/sim/`: 128 simulated invasions |
| `R/05_tables.R` | `output/tables/` |
| `R/06_figures.R` | `output/figures/` |

`R/functions.R` holds the methods (gridding, rarefaction, list-length model, breakpoint search and spread rates). Seeds are set in each script. Steps 02-04 take about 4-5 h on a single core; the slowest runs are those of *Ludwigia grandiflora* (Onagraceae, 2.2 million target-group records). The tables and figures in `output/` are those reported in the paper.

