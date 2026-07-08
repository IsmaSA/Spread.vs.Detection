# =============================================================================
# 02_range_dynamics.R  --  effort-corrected range-dynamics reconstruction
#
# For each species and grid resolution, reconstructs range dynamics under a
# naive model and three independent correction methods:
#   1. naive              cumulative AOO (first-detection cells)
#   2. rarefied           effort-rarefied AOO, continental target-group pool
#   3. rarefied_regional  effort-rarefied AOO, accessible-area-restricted pool
#                         (Barve et al. 2011)
#   4. listlength         list-length detection model (Isaac et al. 2014)
# * project to EPSG:3035, assign integer grid
# * sampled-cell footprint (target-group effort) + relative occupancy
# * invasion-phase breakpoints (segmented) per method
# * radial spread rate (km/yr) per method, with 95% CI
# * cross-method agreement (pairwise Spearman rho; Kendall's W)
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

  # ---- correction method 3: accessible-area-restricted rarefaction ---------
  # same rarefaction machinery, but the effort denominator (target-group pool)
  # is restricted to the species' accessible area (Barve et al. 2011) instead
  # of the whole continent -- tests whether the *scope* of the effort
  # denominator matters, holding the correction *procedure* fixed.
  g_acc <- restrict_target_accessible(f, g, cell_m, cfg$accessible_buffer_km)
  quota_acc <- if (nrow(g_acc) == 0L) quota else
    as.integer(min(quota, stats::median(cumsum(
      rbind(f[, .(year)], g_acc[, .(year)])[, .N, by = year][order(year)]$N))))
  rar_acc <- rarefied_aoo(f, g_acc, years, quota = quota_acc,
                          reps = cfg$rarefy_reps, probs = cfg$rarefy_probs)
  data.table::setnames(rar_acc, c("lwr", "med", "upr"),
                       c("rarr_lwr", "rarr_aoo", "rarr_upr"))
  curve <- merge(curve, rar_acc, by = "year", all.x = TRUE)

  # ---- correction method 4: list-length detection model ---------------------
  # a different correction *paradigm* (detection-probability modelling of
  # presence/absence vs. resampling to constant effort); Isaac et al. (2014).
  ll <- listlength_aoo(f, g, years)
  curve <- merge(curve, ll, by = "year", all.x = TRUE)

  # phase detection per method
  ph_naive <- detect_phases(curve[, .(year, aoo)])
  rar_curve  <- curve[!is.na(rar_aoo),  .(year, aoo = rar_aoo)]
  rarr_curve <- curve[!is.na(rarr_aoo), .(year, aoo = rarr_aoo)]
  ll_curve   <- curve[!is.na(ll_aoo),   .(year, aoo = ll_aoo)]
  ph_rar  <- detect_phases(rar_curve)
  ph_rarr <- detect_phases(rarr_curve)
  ph_ll   <- detect_phases(ll_curve)

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
  w_rarr  <- win(ph_rarr,  rarr_curve)
  w_ll    <- win(ph_ll,    ll_curve)

  sp_naive <- radial_spread(curve[, .(year, aoo)], cell_m, w_naive[1], w_naive[2])
  sp_rar   <- radial_spread(rar_curve,  cell_m, w_rar[1],  w_rar[2])
  sp_rarr  <- radial_spread(rarr_curve, cell_m, w_rarr[1], w_rarr[2])
  sp_ll    <- radial_spread(ll_curve,   cell_m, w_ll[1],   w_ll[2])

  list(
    species = sp, cell_km = cell_km, quota = quota, quota_accessible = quota_acc,
    curve = curve,
    phases = list(naive = ph_naive$breakpoints, rarefied = ph_rar$breakpoints,
                 rarefied_regional = ph_rarr$breakpoints, listlength = ph_ll$breakpoints),
    spread = tibble::tibble(
      species  = sp, cell_km = cell_km,
      method   = c("naive", "rarefied", "rarefied_regional", "listlength"),
      rate_km_yr = c(sp_naive$rate_km_yr, sp_rar$rate_km_yr,
                    sp_rarr$rate_km_yr, sp_ll$rate_km_yr),
      ci_lwr   = c(sp_naive$ci[1], sp_rar$ci[1], sp_rarr$ci[1], sp_ll$ci[1]),
      ci_upr   = c(sp_naive$ci[2], sp_rar$ci[2], sp_rarr$ci[2], sp_ll$ci[2]),
      lag_end  = c(sort(ph_naive$breakpoints)[1], sort(ph_rar$breakpoints)[1],
                  sort(ph_rarr$breakpoints)[1], sort(ph_ll$breakpoints)[1]),
      sat_start= c(sort(ph_naive$breakpoints)[2], sort(ph_rar$breakpoints)[2],
                  sort(ph_rarr$breakpoints)[2], sort(ph_ll$breakpoints)[2])
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

# ---- bias magnitude: naive vs each of the three corrections ---------------
# bias_ratio (naive/rarefied) is the headline metric used throughout the
# Results narrative and Table 3; bias_ratio_regional and bias_ratio_listlength
# repeat the same construction against the other two corrections, so the
# "does naive overestimate spread, and does that scale with citizen-science
# reliance" claim can be checked independently of which correction is used
# as the denominator (not just whether the three corrections agree with each
# other on absolute spread rate, which is the separate §2.6/Fig. 6 question).
# A ratio is only meaningful if its denominator is distinguishable from zero.
# The list-length model in particular can return a fitted rate of ~1e-15
# km/yr (floating-point residue from lm() on an almost perfectly flat
# cumulative-maximum curve, not a genuine near-zero spread rate) for species
# where it detects no year trend; dividing by that produces a "finite" but
# numerically meaningless ratio in the 1e14-1e15 range rather than an error,
# which would silently corrupt any correlation computed on it. Ratios with
# |denominator| below `eps` are set to NA instead.
safe_ratio <- function(num, denom, eps = 0.05) ifelse(abs(denom) < eps, NA_real_, num / denom)

bias_tbl <- spread_all |>
  dplyr::select(species, group, method, rate_km_yr) |>
  tidyr::pivot_wider(names_from = method, values_from = rate_km_yr) |>
  dplyr::mutate(
    bias_ratio            = safe_ratio(naive, rarefied),
    bias_ratio_regional   = safe_ratio(naive, rarefied_regional),
    bias_ratio_listlength = safe_ratio(naive, listlength)
  )

# join mean citizen-science share (post-2010) to test the bias~CS hypothesis
cs <- data.table::fread(file.path(cfg$dir_tab, "citizen_science_share.csv"))
cs_recent <- cs[year >= 2010, .(cs_frac_recent = mean(cs_frac, na.rm = TRUE)),
                by = species]
bias_tbl <- dplyr::left_join(bias_tbl, cs_recent, by = "species")

# correlation test (report r + 95% CI, not just p), one per bias-ratio
# definition; bias_ct is the headline (continental rarefaction) test used
# in the Abstract/Results narrative, the other two are the added checks.
bias_cs_cor <- function(ratio_col) {
  d <- bias_tbl[is.finite(bias_tbl[[ratio_col]]) & is.finite(bias_tbl$cs_frac_recent), ]
  stats::cor.test(d$cs_frac_recent, d[[ratio_col]], method = "pearson")
}
bias_ct            <- bias_cs_cor("bias_ratio")
bias_ct_regional   <- bias_cs_cor("bias_ratio_regional")
bias_ct_listlength <- bias_cs_cor("bias_ratio_listlength")

# ---- leave-one-species-out sensitivity of the bias~citizen-science test ---
# The correlation above uses only 14 points and is visibly influenced by a
# couple of high-leverage species (see Fig. 5 / Results 3.4); this jackknife
# quantifies how much r, its CI, and its P-value shift under each
# single-species exclusion, rather than leaving that as a visual impression.
# Repeated for all three bias-ratio definitions and combined into one long
# table (column `method`) rather than three separate objects.
jackknife_one <- function(ratio_col) {
  d <- bias_tbl[is.finite(bias_tbl[[ratio_col]]) & is.finite(bias_tbl$cs_frac_recent), ]
  purrr::map_dfr(d$species, function(sp_drop) {
    dd <- d[d$species != sp_drop, ]
    ct <- suppressWarnings(stats::cor.test(dd$cs_frac_recent, dd[[ratio_col]], method = "pearson"))
    tibble::tibble(species_dropped = sp_drop,
                   r       = unname(ct$estimate),
                   ci_lwr  = ct$conf.int[1],
                   ci_upr  = ct$conf.int[2],
                   p_value = ct$p.value)
  })
}
jackknife_bias_cor <- dplyr::bind_rows(
  dplyr::mutate(jackknife_one("bias_ratio"),            method = "rarefied"),
  dplyr::mutate(jackknife_one("bias_ratio_regional"),   method = "rarefied_regional"),
  dplyr::mutate(jackknife_one("bias_ratio_listlength"), method = "listlength")
)

# ---- cross-method agreement (robustness to correction method) -------------
# Pairwise Spearman correlations show how similarly any two methods rank
# species by spread rate; Kendall's W summarizes overall concordance among
# the three *corrected* methods (rarefied / rarefied_regional / listlength;
# naive is the uncorrected baseline being tested, not a candidate correction)
# -- i.e. whether conclusions about relative spread rate are robust to the
# choice of correction procedure and effort-denominator scope.
method_mat <- as.matrix(bias_tbl[, c("naive", "rarefied", "rarefied_regional", "listlength")])
rownames(method_mat) <- bias_tbl$species

pairwise_cor <- function(mat) {
  utils::combn(colnames(mat), 2, function(pair) {
    ok <- stats::complete.cases(mat[, pair])
    ct <- suppressWarnings(stats::cor.test(mat[ok, pair[1]], mat[ok, pair[2]],
                                           method = "spearman"))
    tibble::tibble(method_a = pair[1], method_b = pair[2],
                  rho = unname(ct$estimate), p_value = ct$p.value, n = sum(ok))
  }, simplify = FALSE) |> dplyr::bind_rows()
}
method_agreement <- pairwise_cor(method_mat)

corr_mat <- method_mat[, c("rarefied", "rarefied_regional", "listlength")]
corr_mat <- corr_mat[stats::complete.cases(corr_mat), , drop = FALSE]
ranks <- apply(corr_mat, 2, rank, ties.method = "average")
n_kw <- nrow(ranks); m_kw <- ncol(ranks)
Rj <- rowSums(ranks)
S  <- sum((Rj - mean(Rj))^2)
kendall_w <- if (n_kw >= 2L) (12 * S) / (m_kw^2 * (n_kw^3 - n_kw)) else NA_real_

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
readr::write_csv(spread_all,         file.path(cfg$dir_tab, "spread_rates.csv"))
readr::write_csv(bias_tbl,           file.path(cfg$dir_tab, "bias_summary.csv"))
readr::write_csv(sens_all,           file.path(cfg$dir_tab, "sensitivity_grid.csv"))
readr::write_csv(method_agreement,   file.path(cfg$dir_tab, "method_agreement.csv"))
readr::write_csv(jackknife_bias_cor, file.path(cfg$dir_tab, "jackknife_bias_cor.csv"))
data.table::fwrite(curves_all, file.path(cfg$dir_tab, "range_curves.csv"))
saveRDS(list(res = res_primary, curves = curves_all, spread = spread_all,
             bias = bias_tbl, bias_test = bias_ct,
             bias_test_regional = bias_ct_regional,
             bias_test_listlength = bias_ct_listlength,
             sens = sens_all,
             method_agreement = method_agreement, kendall_w = kendall_w,
             jackknife_bias_cor = jackknife_bias_cor),
        file.path(cfg$dir_data, "range_dynamics.rds"))

message("Range dynamics complete | species x resolutions analysed | ",
        "bias r = ", round(unname(bias_ct$estimate), 2),
        " | Kendall's W (corrected methods) = ", round(kendall_w, 2))
