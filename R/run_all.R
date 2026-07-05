# =============================================================================
# run_all.R  --  reproduce the full analysis end to end
#
# Prerequisites:
#   * R 4.5.0, packages listed in 00_setup.R
#   * GBIF credentials in ~/.Renviron (GBIF_USER, GBIF_PWD, GBIF_EMAIL)
#   * run from the project root:  Rscript R/run_all.R
#
# Runtime is dominated by step 01 (GBIF prepares + serves the downloads
# asynchronously; target-group downloads for e.g. Actinopterygii can be large).
# Steps 02-03 run in minutes on the cleaned data. Downloads are cached, so
# re-running only repeats the analysis.
#  =============================================================================

t0 <- Sys.time()
source("R/00_setup.R")                       # config + helpers
source("R/01_data_acquisition.R")            # GBIF download + cleaning  [internet]
source("R/02_range_dynamics.R")              # effort-corrected reconstruction
source("R/03_synthesis_figures_tables.R")    # figures, tables, results_values.json

message("\nDONE in ", round(difftime(Sys.time(), t0, units = "mins"), 1),
        " min. Outputs in output/. Populate manuscript from ",
        "output/results_values.json and output/tables/.")
