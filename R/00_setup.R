# =============================================================================
# 00_setup.R  --  configuration, packages, and helper functions
#
# Project: Separating spread from detection -- effort-corrected reconstruction
#          of freshwater invasion dynamics across Europe from GBIF data.
#
# Author:  I. <lastname>, CSIC
# R:        4.5.0
# CRS:      analysis in ETRS89-LAEA Europe (EPSG:3035, equal-area, metres)
# Style:    native pipe |>, explicit namespacing, colourblind-safe palettes
# =============================================================================

set.seed(20260702)

# ---- packages (pinned versions used during development) ---------------------
# install.packages(c("rgbif","CoordinateCleaner","sf","data.table","dplyr",
#   "tidyr","purrr","stringr","segmented","viridisLite","ggplot2","ggspatial",
#   "rnaturalearth","rnaturalearthdata","readr","tibble","patchwork","scales",
#   "jsonlite"))
suppressPackageStartupMessages({
  library(rgbif)             # >= 3.8.0   GBIF download API
  library(CoordinateCleaner) # >= 3.0.1   coordinate QC
  library(sf)                # >= 1.0-16  vector spatial
  library(data.table)        # >= 1.15.4  memory-efficient tables
  library(dplyr)             # >= 1.1.4
  library(tidyr)             # >= 1.3.1
  library(purrr)             # >= 1.0.2
  library(stringr)           # >= 1.5.1
  library(tibble)            # >= 3.2.1
  library(readr)             # >= 2.1.5
  library(segmented)         # >= 2.1-0   breakpoint regression
  library(ggplot2)           # >= 3.5.1
  library(ggspatial)         # >= 1.1.9   scale bar / north arrow
  library(viridisLite)       # >= 0.4.2
  library(rnaturalearth)     # >= 1.0.1   basemap
  library(patchwork)         # >= 1.2.0   figure composition
  library(scales)            # >= 1.3.0
  library(jsonlite)          # >= 1.8.8   results export
})

# ---- paths ------------------------------------------------------------------
cfg <- list(
  dir_data     = "data",
  dir_dl       = "data/gbif_downloads",   # raw GBIF zips + import cache
  dir_out      = "output",
  dir_fig      = "output/figures",
  dir_tab      = "output/tables",

  # ---- analysis parameters -------------------------------------------------
  crs_analysis = 3035L,                   # ETRS89-LAEA Europe (equal-area, m)
  crs_geo      = 4326L,                   # WGS84 (GBIF native)
  cell_km      = 25,                      # primary grid resolution
  cell_km_sens = c(10, 50),               # sensitivity resolutions
  year_min     = 1950L,
  year_max     = 2024L,
  bin_years    = 1L,                      # temporal bin (annual)

  # Coarse bounding box for the GBIF download predicate (WGS84 degrees). This
  # rectangle is deliberately loose -- it also covers Northern Africa and
  # Turkiye/Asia Minor -- because GBIF download predicates only support
  # simple geometry cheaply; the precise continental-Europe boundary is
  # enforced afterwards in cleaning (see `europe_boundary()` / `restrict_to_europe()`).
  bbox         = c(lon_min = -25, lon_max = 45, lat_min = 34, lat_max = 72),

  # coordinate-uncertainty ceiling (m); records above this are dropped
  coord_unc_max = 10000,

  # effort rarefaction
  rarefy_reps  = 200L,                    # bootstrap replicates
  rarefy_probs = c(0.025, 0.5, 0.975),    # CI + median

  # accessible-area restriction for the regional rarefaction correction
  # (sensu Barve et al. 2011): radius around ever-occupied focal cells used
  # to define the species-specific effort denominator, instead of the whole
  # continent. Km, not a hard biological limit -- checked in sensitivity.
  accessible_buffer_km = 150,

  # basis-of-record types to keep (identification quality proxy)
  keep_basis   = c("HUMAN_OBSERVATION", "MACHINE_OBSERVATION",
                   "PRESERVED_SPECIMEN", "OBSERVATION", "OCCURRENCE"),

  # citizen-science flag: publishingOrg / datasetKey of iNaturalist research grade
  inat_dataset = "50c9509d-22c7-4a22-a47d-8c48425ef4a7"  # iNaturalist RG dataset key
)

dir.create(cfg$dir_dl,  recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$dir_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$dir_tab, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Study species and their target-group effort proxies.
#
# `target_taxon` / `target_rank` define the "target-group background" (sensu
# Phillips et al. 2009) used as the sampling-effort denominator: taxa recorded
# by the same observer community with comparable methods. We use the smallest
# containing higher taxon that is (i) recorded by the same freshwater-sampling
# community and (ii) resolvable in the current GBIF backbone taxonomy -- class
# (Bivalvia), order, or family. Two groups are proxied below the rank that
# would be taxonomically ideal, because the ideal name/rank does not resolve
# via rgbif::name_backbone() (matchType = NONE) despite being valid taxonomy:
#   - Freshwater crayfish: family level (Astacidae / Cambaridae) instead of
#     superfamily Astacoidea, which excludes marine clawed lobsters but exists
#     only in constituent checklists, not the merged backbone.
#   - Ray-finned fish: order level (e.g. Perciformes, Cypriniformes) instead of
#     class Actinopterygii. That class node (usageKey 204) was deleted from the
#     backbone in Nov 2022; current fish species resolve straight from phylum
#     Chordata to order with no class in between, so no class-level proxy is
#     available for this group any more.
# These are transparent, reproducible proxies for freshwater sampling effort,
# not the taxonomically tightest possible grouping.
#
# `target_rank` is ADVISORY metadata: it documents the intended rank and drives
# a post-hoc sanity check in gbif_key(); it is NOT sent as a hard constraint to
# the GBIF matcher (the backbone indexes intermediate ranks inconsistently).
#
# native_region is metadata (literature) for interpretation only; it is NOT used
# in any analysis -- the analysis is GBIF-only.
# =============================================================================
species_tbl <- tibble::tribble(
  ~species,                        ~group,       ~target_taxon,     ~target_rank, ~native_region,
  "Dreissena polymorpha",          "Mollusc",    "Bivalvia",        "class",      "Ponto-Caspian",
  "Dreissena bugensis",            "Mollusc",    "Bivalvia",        "class",      "Ponto-Caspian",
  "Corbicula fluminea",            "Mollusc",    "Bivalvia",        "class",      "E/SE Asia",
  "Procambarus clarkii",           "Crustacean", "Cambaridae",      "family",     "N/C America",
  "Pacifastacus leniusculus",      "Crustacean", "Astacidae",       "family",     "W North America",
  "Faxonius limosus",              "Crustacean", "Cambaridae",      "family",     "E North America",
  "Dikerogammarus villosus",       "Crustacean", "Amphipoda",       "order",      "Ponto-Caspian",
  "Neogobius melanostomus",        "Fish",       "Perciformes",     "order",      "Ponto-Caspian",
  "Pseudorasbora parva",           "Fish",       "Cypriniformes",   "order",      "E Asia",
  "Lepomis gibbosus",              "Fish",       "Perciformes",     "order",      "E North America",
  "Gambusia holbrooki",            "Fish",       "Cyprinodontiformes","order",    "SE North America",
  "Ameiurus melas",                "Fish",       "Siluriformes",    "order",      "E North America",
  "Elodea nuttallii",              "Plant",      "Hydrocharitaceae","family",     "North America",
  "Ludwigia grandiflora",          "Plant",      "Onagraceae",      "family",     "S America"
)

# =============================================================================
# Helper functions
# =============================================================================

# null/NA-coalescing helper (base R >= 4.4 exports `%||%`; defined here for
# portability to older R and to also treat length-0 / all-NA as "empty").
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L || all(is.na(a))) b else a

#' Resolve a taxon name to its accepted GBIF backbone usageKey.
#'
#' `rank` is ADVISORY: it drives a post-hoc sanity check and the diagnostic
#' message, but is NOT sent as a hard constraint to the matcher. GBIF indexes
#' intermediate ranks (infraorder, superfamily) inconsistently, so a
#' rank-constrained query can return matchType "NONE" even for a valid name.
#' We match on name (+ optional kingdom to break homonyms), then verify.
#'
#' @param name           taxon name at any rank
#' @param rank           expected rank (advisory; sanity check + message only)
#' @param kingdom        optional kingdom hint to disambiguate cross-kingdom homonyms
#' @param min_confidence minimum GBIF match confidence (0-100) before warning
#' @param verbose        print a one-line resolution report per taxon
#' @return integer usageKey (accepted backbone key)
gbif_key <- function(name, rank = "species", kingdom = NULL,
                     min_confidence = 90L, verbose = TRUE) {
  args <- list(name = name, strict = FALSE)
  if (!is.null(kingdom)) args$kingdom <- kingdom
  m <- do.call(rgbif::name_backbone, args)

  key <- suppressWarnings(as.integer(m$usageKey %||% NA_integer_))
  if (isTRUE(m$matchType == "NONE") || is.na(key))
    stop("No GBIF backbone match for '", name, "' (matchType = ",
         as.character(m$matchType %||% "NA"),
         "). Check spelling, drop the rank hint, or try a neighbouring rank.")

  matched_rank <- tolower(as.character(m$rank %||% NA_character_))
  conf         <- suppressWarnings(as.integer(m$confidence %||% NA_integer_))

  if (!is.na(conf) && conf < min_confidence)
    warning("Low-confidence GBIF match for '", name, "': confidence ", conf,
            " (< ", min_confidence, ").")
  want_species <- identical(tolower(rank), "species")
  if (want_species && !is.na(matched_rank) && matched_rank != "species")
    warning("Expected a species for '", name, "' but GBIF returned rank '",
            matched_rank, "'.")
  if (!want_species && !is.na(matched_rank) && matched_rank %in% c("species", "genus"))
    warning("Expected a higher taxon for '", name, "' but GBIF returned '",
            matched_rank, "' (key ", key, "); the name may be unmatched at its true rank.")

  if (verbose)
    message(sprintf("  GBIF: %-20s key=%-9s rank=%-11s conf=%s",
                    name, key, matched_rank, as.character(conf)))
  key
}

#' Precise continental-Europe boundary (dissolved country polygons), used to
#' drop records that fall inside the coarse download bounding box (`cfg$bbox`)
#' but outside the study region -- chiefly Northern Africa and Turkiye/Asia
#' Minor. Natural Earth's continent classification puts Turkiye in Asia and
#' Morocco/Algeria/Tunisia/Libya/Egypt in Africa, so filtering to
#' continent == "Europe" excludes both without a manual country-code list.
#' @return a single (multi)polygon sfc, WGS84 (EPSG:4326)
europe_boundary <- function() {
  rnaturalearth::ne_countries(scale = "medium", continent = "Europe",
                              returnclass = "sf") |>
    sf::st_make_valid() |>
    sf::st_union()
}

#' Restrict a cleaned occurrence data.table (lon/lat in WGS84) to points
#' falling within `europe_poly` (see `europe_boundary()`).
#' @param dt data.table with columns lon, lat
#' @param europe_poly sfc polygon, WGS84 (EPSG:4326)
restrict_to_europe <- function(dt, europe_poly) {
  if (nrow(dt) == 0L) return(dt)
  pts <- sf::st_as_sf(dt, coords = c("lon", "lat"), crs = cfg$crs_geo, remove = FALSE)
  inside <- lengths(sf::st_intersects(pts, europe_poly)) > 0
  dt[inside]
}

#' Assign projected points to an integer grid (memory-efficient, no grid object).
#' @param dt data.table with X, Y in metres (EPSG:3035)
#' @param cell_m grid cell size in metres
#' @return dt with integer columns gx, gy and a character cell id
assign_grid <- function(dt, cell_m) {
  dt <- data.table::copy(dt)
  dt[, gx := as.integer(floor(X / cell_m))]
  dt[, gy := as.integer(floor(Y / cell_m))]
  dt[, cell := paste(gx, gy, sep = "_")]
  dt[]
}

#' Naive cumulative area of occupancy: distinct cells occupied by first-detection
#' year. Returns one row per year from year_min..year_max.
naive_aoo <- function(dt_focal, years) {
  first_seen <- dt_focal[, .(first_year = min(year)), by = cell]
  n_new <- first_seen[, .N, by = first_year][order(first_year)]
  out <- data.table::data.table(year = years)
  out <- merge(out, n_new, by.x = "year", by.y = "first_year", all.x = TRUE)
  out[is.na(N), N := 0L]
  out[, aoo := cumsum(N)]
  out[, .(year, aoo)]
}

#' Cumulative number of *sampled* cells (target-group effort footprint).
sampled_cells <- function(dt_target, years) {
  first_samp <- dt_target[, .(first_year = min(year)), by = cell]
  n_new <- first_samp[, .N, by = first_year][order(first_year)]
  out <- data.table::data.table(year = years)
  out <- merge(out, n_new, by.x = "year", by.y = "first_year", all.x = TRUE)
  out[is.na(N), N := 0L]
  out[, sampled := cumsum(N)]
  out[, .(year, sampled)]
}

#' Effort-rarefied AOO. For each evaluation year the cumulative record pool
#' (focal + target) is subsampled to a fixed size `quota`; distinct focal cells
#' in the subsample are counted. Holding sample size constant removes the shared
#' effort trend, so change over years reflects true occupancy (cf. coverage-based
#' rarefaction). Returns median and CI across `reps`.
rarefied_aoo <- function(dt_focal, dt_target, years, quota,
                         reps = 200L, probs = c(0.025, 0.5, 0.975)) {
  # pooled records: is_focal flag + cell + year
  pool <- data.table::rbindlist(list(
    dt_focal[,  .(cell, year, is_focal = TRUE)],
    dt_target[, .(cell, year, is_focal = FALSE)]
  ))
  data.table::setkey(pool, year)
  eval_years <- years[years >= min(dt_focal$year)]

  res <- vector("list", length(eval_years))
  for (i in seq_along(eval_years)) {
    y  <- eval_years[i]
    cp <- pool[year <= y]
    n  <- nrow(cp)
    if (n < quota) { res[[i]] <- data.table::data.table(year = y,
                       lwr = NA_real_, med = NA_real_, upr = NA_real_); next }
    vals <- vapply(seq_len(reps), function(b) {
      idx <- sample.int(n, quota, replace = FALSE)
      data.table::uniqueN(cp[idx][is_focal == TRUE, cell])
    }, numeric(1))
    q <- stats::quantile(vals, probs = probs, names = FALSE)
    res[[i]] <- data.table::data.table(year = y, lwr = q[1], med = q[2], upr = q[3])
  }
  data.table::rbindlist(res)
}

#' Restrict a target-group point set to the "accessible area" around a focal
#' species' realized range (sensu Barve et al. 2011): target-group cells whose
#' centroid lies within `buffer_km` of ANY ever-occupied focal cell. Used as
#' an alternative, species-specific effort denominator (vs. the whole
#' continental target group) for the regional rarefaction correction.
#' @param dt_focal  focal-species gridded points (has gx, gy from assign_grid)
#' @param dt_target target-group gridded points (has gx, gy from assign_grid)
#' @param cell_m    grid cell size in metres (must match assign_grid call)
#' @param buffer_km buffer radius, km
#' @return dt_target filtered to rows within the accessible area
restrict_target_accessible <- function(dt_focal, dt_target, cell_m, buffer_km) {
  occ <- unique(dt_focal[, .(gx, gy)])
  if (nrow(occ) == 0L) return(dt_target[0])
  occ_cx <- (occ$gx + 0.5) * cell_m
  occ_cy <- (occ$gy + 0.5) * cell_m

  tgt_cells <- unique(dt_target[, .(cell, gx, gy)])
  tgt_cells[, `:=`(cx = (gx + 0.5) * cell_m, cy = (gy + 0.5) * cell_m)]
  buffer_m <- buffer_km * 1000

  # min distance from each target cell to the nearest ever-occupied focal cell
  min_dist <- vapply(seq_len(nrow(tgt_cells)), function(i) {
    min(sqrt((occ_cx - tgt_cells$cx[i])^2 + (occ_cy - tgt_cells$cy[i])^2))
  }, numeric(1))
  keep_cells <- tgt_cells$cell[min_dist <= buffer_m]
  dt_target[cell %in% keep_cells]
}

#' Detection-corrected AOO via a list-length model (Szabo et al. 2010;
#' Isaac et al. 2014). Treats each cell-year with >=1 target-group record as
#' a "checklist" of length = number of target-group records (the recording
#' effort in that cell-year), and models focal-species detection as a
#' function of list length and year:
#'   presence ~ log1p(list_length) + factor(year)
#' Predicting at a fixed reference list length for every year removes the
#' effect of rising effort, giving an occupancy-probability trend that is
#' converted to cell units (probability x cells ever sampled) and forced
#' non-decreasing (cummax), since AOO is a cumulative quantity by construction.
#' @param dt_focal  focal-species gridded points (cell, year)
#' @param dt_target target-group gridded points (cell, year)
#' @param years     full evaluation year sequence
#' @return data.table(year, ll_aoo)
listlength_aoo <- function(dt_focal, dt_target, years) {
  na_out <- data.table::data.table(year = years, ll_aoo = NA_real_)

  tgt_yr <- dt_target[, .(list_length = .N), by = .(cell, year)]
  foc_yr <- unique(dt_focal[, .(cell, year)])[, present := 1L]
  chk <- merge(tgt_yr, foc_yr, by = c("cell", "year"), all.x = TRUE)
  chk[is.na(present), present := 0L]
  chk <- chk[year %in% years]

  if (data.table::uniqueN(chk$year) < 5L || sum(chk$present) < 5L) return(na_out)

  fit <- tryCatch(
    stats::glm(present ~ log1p(list_length) + factor(year), data = chk,
              family = stats::binomial()),
    error = function(e) NULL)
  if (is.null(fit)) return(na_out)

  ref_ll  <- stats::median(chk$list_length)
  eval_yr <- sort(unique(chk$year))
  pred <- data.table::data.table(year = eval_yr, list_length = ref_ll)
  pred[, prob := stats::predict(fit, newdata = pred, type = "response")]
  n_cells_total <- data.table::uniqueN(chk$cell)
  pred[, ll_aoo := cummax(prob * n_cells_total)]

  out <- merge(data.table::data.table(year = years), pred[, .(year, ll_aoo)],
              by = "year", all.x = TRUE)
  out[, ll_aoo := data.table::nafill(ll_aoo, type = "locf")]
  out[]
}

#' Radial spread rate (km/yr) via the diffusion approximation:
#' sqrt(occupied area / pi) regressed on year over the expansion phase.
#' @param aoo_dt data.table(year, aoo) with aoo = number of occupied cells
#' @param cell_m cell size (m); area per cell = cell_m^2
#' @param yr_start,yr_end expansion-phase bounds
radial_spread <- function(aoo_dt, cell_m, yr_start, yr_end) {
  d <- aoo_dt[year >= yr_start & year <= yr_end & aoo > 0]
  if (nrow(d) < 3) return(list(rate_km_yr = NA_real_, ci = c(NA, NA), n = nrow(d)))
  d[, area_km2 := aoo * (cell_m / 1000)^2]
  d[, radius_km := sqrt(area_km2 / pi)]
  fit <- stats::lm(radius_km ~ year, data = d)
  ci  <- stats::confint(fit, "year", level = 0.95)
  list(rate_km_yr = unname(coef(fit)["year"]),
       ci = as.numeric(ci), n = nrow(d), model = fit)
}

#' Detect invasion phases (lag / expansion / saturation) as breakpoints in
#' log(cumulative AOO + 1) vs year using segmented regression (2 breakpoints).
detect_phases <- function(aoo_dt, npsi = 2) {
  d <- aoo_dt[aoo > 0]
  if (nrow(d) < (npsi + 3)) return(list(breakpoints = rep(NA_real_, npsi)))
  d[, y := log1p(aoo)]
  lm0 <- stats::lm(y ~ year, data = d)
  sg  <- tryCatch(
    segmented::segmented(lm0, seg.Z = ~year, npsi = npsi),
    error = function(e) NULL)
  if (is.null(sg) || is.null(sg$psi))
    return(list(breakpoints = rep(NA_real_, npsi), model = NULL))
  list(breakpoints = round(sg$psi[, "Est."]),
       slopes = segmented::slope(sg)$year[, "Est."], model = sg)
}

#' Gompertz fit to a cumulative curve -> asymptote (K) and inflection year.
fit_gompertz <- function(aoo_dt) {
  d <- aoo_dt[aoo > 0]
  if (nrow(d) < 5) return(NULL)
  st <- list(K = max(d$aoo) * 1.1, xmid = stats::median(d$year), scal = 5)
  tryCatch(
    stats::nls(aoo ~ K * exp(-exp(-(year - xmid) / scal)),
               data = d, start = st,
               control = stats::nls.control(maxiter = 500, warnOnly = TRUE)),
    error = function(e) NULL)
}

# colourblind-safe palette for the taxonomic groups
group_pal <- c(Mollusc = "#009E73", Crustacean = "#E69F00",
               Fish = "#0072B2", Plant = "#CC79A7")

# colourblind-safe palette for the four spread-rate correction methods
method_pal <- c(naive = "#D55E00", rarefied = "#0072B2",
                rarefied_regional = "#009E73", listlength = "#CC79A7")
method_lab <- c(naive = "Naive", rarefied = "Rarefied (continental)",
               rarefied_regional = "Rarefied (accessible area)",
               listlength = "List-length model")

theme_inv <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(face = "italic"),
                   legend.position = "bottom")
}

message("Setup complete | ", nrow(species_tbl), " species | grid ",
        cfg$cell_km, " km | CRS EPSG:", cfg$crs_analysis)
