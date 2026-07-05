# =============================================================================
# 02_range_dynamics.R  --  effort-corrected range-dynamics reconstruction
#
# For each species and grid resolution:
#   * project to EPSG:3035, assign integer grid
#   * naive cumulative AOO (first-detection cells)
#   * sampled-cell footprint (target-group effort)
#   * relative occupancy = AOO / sampled
#   * effort-rarefied AOO (constant sample size)
#   * invasion-phase breakpoints (segmented) on naive vs rarefied curves
#   * radial spread rate (km/yr), naive vs rarefied, with 95% CI
# =============================================================================

source("R/00_setup.R")

occ <- readRDS(file.path(cfg$dir_data, "occ_clean.rds"))
species_tbl <- occ$species
years <- cfg$year_min:cfg$year_max

# project WGS84 -> EPSG:3035, return data.table(cell coords in metres)
project_dt <- function(dt) {
  p <- sf::st_as_sf(dt, coords = c("lon", "lat"), crs = cfg$crs_geo) |>
    sf::st_transform(cfg$crs_analysis)
  xy <- sf::st_coordinates(p)
  data.table::data.table(X = xy[, 1], Y = xy[, 2], year = dt$year)
}

focal_xy  <- purrr::map(occ$focal,  project_dt)
target_xy <- purrr::map(occ$target, project_dt)

# map each species to its target-group projected points
sp_target <- setNames(species_tbl$target_taxon, species_tbl$species)

# -----------------------------------------------------------------------------
# analyse one species at one grid resolution
# -----------------------------------------------------------------------------
analyse_species <- function(sp, cell_km) {
  cell_m <- cell_km * 1000
  f <- assign_grid(focal_xy[[sp]],  cell_m)
  g <- assign_grid(target_xy[[sp_target[[sp]]]], cell_m)

  naive <- naive_aoo(f, years)
  samp  <- sampled_cells(g, years)
  curve <- merge(naive, samp, by = "year")
  curve[, rel_occ := ifelse(sampled > 0, aoo / sampled, NA_real_)]

  # rarefaction quota: median annual cumulative pool size, capped for speed
  pool_by_year <- rbind(f[, .(year)], g[, .(year)])[, .N, by = year][order(year)]
  quota <- as.integer(min(50000, stats::median(cumsum(pool_by_year$N))))
  rar <- rarefied_aoo(f, g, years, quota = quota,
                      reps = cfg$rarefy_reps, probs = cfg$rarefy_probs)
  data.table::setnames(rar, c("lwr", "med", "upr"),
                       c("rar_lwr", "rar_aoo", "rar_upr"))
  curve <- merge(curve, rar, by = "year", all.x = TRUE)

  # phase detection on naive and rarefied curves
  ph_naive <- detect_phases(curve[, .(year, aoo)])
  rar_curve <- curve[!is.na(rar_aoo), .(year, aoo = rar_aoo)]
  ph_rar    <- detect_phases(rar_curve)

  # expansion window = between the two breakpoints (fallback to central 60%)
  win <- function(ph, dd) {
    bp <- ph$breakpoints
    if (all(is.na(bp))) {
      yr <- range(dd[aoo > 0, year]); span <- diff(yr)
      return(c(yr[1] + 0.2 * span, yr[2] - 0.2 * span))
    }
    sort(bp)[1:2]
  }
  w_naive <- win(ph_naive, curve[, .(year, aoo)])
  w_rar   <- win(ph_rar,   rar_curve)

  sp_naive <- radial_spread(curve[, .(year, aoo)], cell_m, w_naive[1], w_naive[2])
  sp_rar   <- radial_spread(rar_curve, cell_m, w_rar[1], w_rar[2])

  list(
    species = sp, cell_km = cell_km, quota = quota, curve = curve,
    phases = list(naive = ph_naive$breakpoints, rarefied = ph_rar$breakpoints),
    spread = tibble::tibble(
      species  = sp, cell_km = cell_km,
      method   = c("naive", "rarefied"),
      rate_km_yr = c(sp_naive$rate_km_yr, sp_rar$rate_km_yr),
      ci_lwr   = c(sp_naive$ci[1], sp_rar$ci[1]),
      ci_upr   = c(sp_naive$ci[2], sp_rar$ci[2]),
      lag_end  = c(sort(ph_naive$breakpoints)[1], sort(ph_rar$breakpoints)[1]),
      sat_start= c(sort(ph_naive$breakpoints)[2], sort(ph_rar$breakpoints)[2])
    )
  )
}

# ---- run at primary resolution ---------------------------------------------
res_primary <- purrr::map(species_tbl$species, analyse_species,
                          cell_km = cfg$cell_km) |>
  setNames(species_tbl$species)

curves_all <- purrr::imap(res_primary, ~ .x$curve[, species := .y]) |>
  data.table::rbindlist(fill = TRUE)
spread_all <- purrr::map(res_primary, "spread") |> dplyr::bind_rows() |>
  dplyr::left_join(dplyr::select(species_tbl, species, group, native_region),
                   by = "species")

# ---- bias magnitude: naive vs rarefied recent expansion --------------------
bias_tbl <- spread_all |>
  dplyr::select(species, group, method, rate_km_yr) |>
  tidyr::pivot_wider(names_from = method, values_from = rate_km_yr) |>
  dplyr::mutate(bias_ratio = naive / rarefied)

# join mean citizen-science share (post-2010) to test the bias~CS hypothesis
cs <- data.table::fread(file.path(cfg$dir_tab, "citizen_science_share.csv"))
cs_recent <- cs[year >= 2010, .(cs_frac_recent = mean(cs_frac, na.rm = TRUE)),
                by = species]
bias_tbl <- dplyr::left_join(bias_tbl, cs_recent, by = "species")

# correlation test (report r + 95% CI, not just p)
bias_ct <- with(bias_tbl[is.finite(bias_tbl$bias_ratio), ],
                stats::cor.test(cs_frac_recent, bias_ratio, method = "pearson"))

# ---- sensitivity across grid resolutions -----------------------------------
res_sens <- purrr::map(cfg$cell_km_sens, function(ck) {
  purrr::map(species_tbl$species, analyse_species, cell_km = ck) |>
    purrr::map("spread") |> dplyr::bind_rows()
}) |> dplyr::bind_rows()
sens_all <- dplyr::bind_rows(
  dplyr::mutate(dplyr::select(spread_all, species, cell_km, method, rate_km_yr), set = "primary"),
  dplyr::mutate(res_sens[, c("species", "cell_km", "method", "rate_km_yr")], set = "sensitivity")
)

# ---- persist ---------------------------------------------------------------
readr::write_csv(spread_all, file.path(cfg$dir_tab, "spread_rates.csv"))
readr::write_csv(bias_tbl,   file.path(cfg$dir_tab, "bias_summary.csv"))
readr::write_csv(sens_all,   file.path(cfg$dir_tab, "sensitivity_grid.csv"))
data.table::fwrite(curves_all, file.path(cfg$dir_tab, "range_curves.csv"))
saveRDS(list(res = res_primary, curves = curves_all, spread = spread_all,
             bias = bias_tbl, bias_test = bias_ct, sens = sens_all),
        file.path(cfg$dir_data, "range_dynamics.rds"))

message("Range dynamics complete | species x resolutions analysed | ",
        "bias r = ", round(unname(bias_ct$estimate), 2))
