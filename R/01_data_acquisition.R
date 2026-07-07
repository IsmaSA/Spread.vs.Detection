# =============================================================================
# 01_data_acquisition.R  --  GBIF downloads + coordinate cleaning
#
# Requires GBIF credentials in ~/.Renviron:
#   GBIF_USER=...  GBIF_PWD=...  GBIF_EMAIL=...
#
# Two download tiers per species:
#   (1) focal    : occurrences of the invader
#   (2) target   : occurrences of the containing class/order (effort proxy)
#
# All downloads are asynchronous and return a citable DOI (recorded for the
# Data Availability statement). Downloads are cached; re-running skips completed
# requests. This is the ONLY step that needs internet + GBIF login.
# =============================================================================
usethis::edit_r_environ()
source("R/00_setup.R")

# ---- 1. resolve backbone keys ----------------------------------------------
species_tbl <- species_tbl |>
  dplyr::mutate(
    focal_key  = purrr::map_int(species, gbif_key),
    target_key = purrr::map2_int(target_taxon, target_rank,
                                 ~ gbif_key(.x, rank = .y))
  )
readr::write_csv(species_tbl, file.path(cfg$dir_tab, "species_keys.csv"))

# ---- 2. shared download predicate (Europe bbox + years + quality) ----------
bb <- cfg$bbox
base_preds <- function(taxon_key) {
  list(
    rgbif::pred("taxonKey", taxon_key),
    rgbif::pred("hasCoordinate", TRUE),
    rgbif::pred("hasGeospatialIssue", FALSE),
    rgbif::pred("occurrenceStatus", "PRESENT"),
    rgbif::pred_gte("year", cfg$year_min),
    rgbif::pred_lte("year", cfg$year_max),
    rgbif::pred_gte("decimalLatitude",  bb[["lat_min"]]),
    rgbif::pred_lte("decimalLatitude",  bb[["lat_max"]]),
    rgbif::pred_gte("decimalLongitude", bb[["lon_min"]]),
    rgbif::pred_lte("decimalLongitude", bb[["lon_max"]])
  )
}

#' Submit (or reuse) a GBIF download and import selected columns.
#' Three cache files per tag: <tag>.csv (imported data), <tag>_dl.rds (submission
#' object, for resuming an in-progress request), <tag>_meta.rds (citable DOI).
get_gbif <- function(taxon_key, tag) {
  csv_file  <- file.path(cfg$dir_dl, paste0(tag, ".csv"))
  dl_file   <- file.path(cfg$dir_dl, paste0(tag, "_dl.rds"))
  meta_file <- file.path(cfg$dir_dl, paste0(tag, "_meta.rds"))
  if (file.exists(csv_file)) {
    message("cache hit: ", tag)
    return(data.table::fread(csv_file))
  }
  dl <- if (file.exists(dl_file)) readRDS(dl_file) else {
    d <- do.call(rgbif::occ_download,
                 c(base_preds(taxon_key), list(format = "SIMPLE_CSV")))
    saveRDS(d, dl_file); d
  }
  rgbif::occ_download_wait(dl)                       # poll until READY
  z <- as.character(rgbif::occ_download_get(dl, path = cfg$dir_dl, overwrite = TRUE))

  # extract with base R (portable: no dependency on an external `unzip` binary,
  # which isn't on PATH by default on Windows)
  exdir <- file.path(tempdir(), tag)
  utils::unzip(z, exdir = exdir)
  csv_in_zip <- list.files(exdir, pattern = "\\.csv$", full.names = TRUE)[1]

  # memory-efficient import: only the columns we need
  keep_cols <- c("gbifID", "species", "decimalLatitude", "decimalLongitude",
                 "coordinateUncertaintyInMeters", "year", "basisOfRecord",
                 "datasetKey", "countryCode")
  dt <- data.table::fread(csv_in_zip, select = keep_cols, quote = "",
                          showProgress = FALSE)
  unlink(exdir, recursive = TRUE)
  data.table::setnames(dt, c("decimalLatitude", "decimalLongitude"),
                       c("lat", "lon"))
  saveRDS(list(doi = rgbif::occ_download_meta(dl)$doi), meta_file)
  data.table::fwrite(dt, csv_file)
  dt
}

# ---- 3. download focal + target for every species --------------------------
focal_raw  <- purrr::pmap(
  list(species_tbl$focal_key, species_tbl$species),
  function(k, sp) get_gbif(k, paste0("focal_", stringr::str_replace_all(sp, " ", "_")))
) |> setNames(species_tbl$species)

# target groups are shared across species with the same proxy -> download once
targets <- species_tbl |> dplyr::distinct(target_taxon, target_key)
target_raw <- purrr::map2(
  targets$target_key, targets$target_taxon,
  ~ get_gbif(.x, paste0("target_", stringr::str_replace_all(.y, " ", "_")))
) |> setNames(targets$target_taxon)

# ---- 4. coordinate cleaning (CoordinateCleaner) ----------------------------
# Freshwater taxa may occur legitimately close to the coast, so cc_sea() is NOT
# applied. We flag capitals, centroids, institutions, GBIF HQ, zeros, equal
# lat/lon, and enforce the coordinate-uncertainty ceiling. The download bbox
# (cfg$bbox) is a coarse rectangle that also covers Northern Africa and
# Turkiye/Asia Minor; restrict_to_europe() enforces the true study region
# (continental Europe) precisely, using a dissolved Natural Earth boundary
# computed once and reused across species.
europe_poly <- europe_boundary()

clean_occ <- function(dt) {
  dt <- dt[!is.na(lat) & !is.na(lon) & !is.na(year)]
  dt <- dt[basisOfRecord %in% cfg$keep_basis]
  dt <- dt[is.na(coordinateUncertaintyInMeters) |
             coordinateUncertaintyInMeters <= cfg$coord_unc_max]
  dt <- unique(dt, by = c("lat", "lon", "year", "species"))
  dt <- restrict_to_europe(dt, europe_poly)
  df <- as.data.frame(dt)
  flags <- CoordinateCleaner::clean_coordinates(
    x = df, lon = "lon", lat = "lat",
    species = if ("species" %in% names(df)) "species" else NULL,
    tests = c("capitals", "centroids", "equal", "gbif", "institutions", "zeros"),
    value = "flagged", verbose = FALSE)
  data.table::as.data.table(df[flags, , drop = FALSE])
}

focal_clean  <- purrr::map(focal_raw,  clean_occ)
target_clean <- purrr::map(target_raw, clean_occ)

# ---- 5. citizen-science share (for the bias analysis) ----------------------
cs_share <- function(dt) {
  dt[, .(n = .N,
         n_inat = sum(datasetKey == cfg$inat_dataset, na.rm = TRUE)),
     by = year][, cs_frac := n_inat / n][order(year)]
}
focal_cs <- purrr::imap(focal_clean, ~ .x |> cs_share() |>
                          dplyr::mutate(species = .y)) |>
  data.table::rbindlist()
readr::write_csv(focal_cs, file.path(cfg$dir_tab, "citizen_science_share.csv"))

# ---- 6. record download DOIs for the manuscript ----------------------------
doi_tbl <- tibble::tibble(
  tag = c(paste0("focal_", species_tbl$species),
          paste0("target_", targets$target_taxon)),
  doi = c(purrr::map_chr(species_tbl$species, ~ {
            m <- file.path(cfg$dir_dl,
                   paste0("focal_", stringr::str_replace_all(.x, " ", "_"), "_meta.rds"))
            if (file.exists(m)) readRDS(m)$doi %||% NA_character_ else NA_character_ }),
          purrr::map_chr(targets$target_taxon, ~ {
            m <- file.path(cfg$dir_dl,
                   paste0("target_", stringr::str_replace_all(.x, " ", "_"), "_meta.rds"))
            if (file.exists(m)) readRDS(m)$doi %||% NA_character_ else NA_character_ }))
)
readr::write_csv(doi_tbl, file.path(cfg$dir_tab, "gbif_download_dois.csv"))

saveRDS(list(focal = focal_clean, target = target_clean, species = species_tbl),
        file.path(cfg$dir_data, "occ_clean.rds"))

message("Data acquisition complete | focal spp: ", length(focal_clean),
        " | target groups: ", length(target_clean))
