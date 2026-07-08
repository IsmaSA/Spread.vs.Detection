# =============================================================================
# 03_synthesis_figures_tables.R  --  publication-ready figures + tables
#
# Figures (300 dpi, colourblind-safe viridis/Okabe-Ito, EPSG:3035 maps):
#   Fig 1  occurrence density of study species across Europe
#   Fig 2  sampling effort through time (the confound) + citizen-science share
#   Fig 3  naive vs effort-corrected range curves per species
#   Fig 4  spread velocity, naive vs three correction methods, with 95% CI
#   Fig 5  cross-species synthesis: bias ratio vs citizen-science share
#   Fig 6  agreement among correction methods (robustness check)
#   Fig S2 leave-one-species-out sensitivity of the bias~CS correlation
# Tables:
#   Table 1  species, group, native region, n records, span, DOI
#   Table 2  invasion-phase breakpoints (naive vs three correction methods)
#   Table 3  spread rates (km/yr) with CIs, all four methods
# Plus results_values.json  (numbers the manuscript template pulls in).
# =============================================================================

source("R/00_setup.R")

occ  <- readRDS(file.path(cfg$dir_data, "occ_clean.rds"))
rd   <- readRDS(file.path(cfg$dir_data, "range_dynamics.rds"))
sp_t <- occ$species
curves <- data.table::as.data.table(rd$curves)
ggsave2 <- function(f, p, w, h) ggplot2::ggsave(file.path(cfg$dir_fig, f), p,
                                                width = w, height = h, dpi = 300)

# =============================================================================
# Figure 1 -- occupied grid cells across Europe
# =============================================================================
map_bbox <- sf::st_bbox(c(xmin = 2.5e6, xmax = 6e6, ymin = 1.5e6, ymax = 5.4e6),
                        crs = cfg$crs_analysis)

europe <- rnaturalearth::ne_countries(scale = "medium", continent = "Europe",
                                       returnclass = "sf") |>
  sf::st_transform(cfg$crs_analysis)

# Occupied cells on the *same* grid used for the AOO analysis (cfg$cell_km),
# rather than raw occurrence points -- this ties the figure directly to what
# "area of occupancy" means in this paper, and plots at most a few thousand
# cells instead of ~253k individual points (which previously needed
# rasterising just to keep the SVG file size and element count manageable).
cell_m  <- cfg$cell_km * 1000
grid_sf <- purrr::imap(occ$focal, ~ .x[, .(lon, lat)][, species := .y]) |>
  data.table::rbindlist() |>
  sf::st_as_sf(coords = c("lon", "lat"), crs = cfg$crs_geo) |>
  sf::st_transform(cfg$crs_analysis)
xy <- sf::st_coordinates(grid_sf)
grid_dt <- data.table::data.table(species = grid_sf$species, X = xy[, 1], Y = xy[, 2]) |>
  assign_grid(cell_m)
grid_dt <- dplyr::left_join(grid_dt, dplyr::select(sp_t, species, group), by = "species")
grid_dt[, `:=`(cx = (gx + 0.5) * cell_m, cy = (gy + 0.5) * cell_m)]

# panel (a): cells ever occupied by ANY species in a group, pooled across
# species. Row order is randomized so overlapping group tiles (a cell
# invaded by more than one group) don't systematically favour whichever
# group happens to be drawn last.
occ_group <- unique(grid_dt[, .(group, cx, cy)])
occ_group <- occ_group[sample.int(nrow(occ_group)), ]

# panel (b): cells occupied per species (no overlap within a facet panel).
occ_species <- unique(grid_dt[, .(species, group, cx, cy)])

# ---- panel (a): all species pooled, occupied cells by taxonomic group -----
fig1a <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = europe, fill = "grey96", colour = "grey80", linewidth = 0.2) +
  ggplot2::geom_tile(data = occ_group, ggplot2::aes(cx, cy, fill = group),
                     width = cell_m, height = cell_m, alpha = 0.6) +
  ggplot2::scale_fill_manual(values = group_pal, name = "Group") +
  ggplot2::coord_sf(crs = sf::st_crs(cfg$crs_analysis),
                    xlim = map_bbox[c("xmin", "xmax")], ylim = map_bbox[c("ymin", "ymax")]) +
  ggspatial::annotation_scale(location = "bl", height = grid::unit(0.15, "cm")) +
  theme_inv(10) +
  ggplot2::labs(subtitle = "(a) All 14 species pooled, occupied cells by taxonomic group")

# ---- panel (b): occupied cells faceted by species --------------------------
fig1b <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = europe, fill = "grey76", colour = "grey80", linewidth = 0.2) +
  ggplot2::geom_tile(data = occ_species, ggplot2::aes(cx, cy, fill = group),
                     width = cell_m, height = cell_m, alpha = 0.85) +
  ggplot2::scale_fill_manual(values = group_pal, name = "Group") +
  ggplot2::coord_sf(crs = sf::st_crs(cfg$crs_analysis),
                    xlim = map_bbox[c("xmin", "xmax")], ylim = map_bbox[c("ymin", "ymax")]) +
  ggplot2::facet_wrap(~species, ncol = 4) +
  theme_inv(9) +
  ggplot2::theme(strip.text = ggplot2::element_text(size = 14, face = "italic"))

fig1 <- patchwork::wrap_plots(fig1a, fig1b, ncol = 1, heights = c(0.3, 0.7)) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title = paste0("Cells occupied by study invaders (", cfg$cell_km, " km grid; 1950-2024)"),
    caption = "ETRS89-LAEA (EPSG:3035). A cell is shaded if the species/group was recorded there at least once.")

ggsave2("Fig1_occurrence_maps.svg", fig1b, 11, 13)

# =============================================================================
# Figure 2 -- sampling effort through time + citizen-science share
# =============================================================================
effort_ts <- curves[, .(sampled = max(sampled, na.rm = TRUE)), by = .(species, year)]
cs <- data.table::fread(file.path(cfg$dir_tab, "citizen_science_share.csv"))
x_breaks <- seq(1950, 2025, by = 20)
xx_breaks <- seq(1990, 2025, by = 10)

p2a <- ggplot2::ggplot(effort_ts, ggplot2::aes(year, sampled, colour = species)) +
  ggplot2::geom_line(alpha = 0.85, linewidth = 0.6) +
  ggplot2::scale_y_log10() +
  ggplot2::scale_colour_viridis_d(name = "Species") +
  theme_inv() + theme_bw() +  ggplot2::scale_x_continuous(
    limits = c(1950, 2025),
    breaks = x_breaks
  ) + ggplot2::theme(
    axis.title = ggplot2::element_text(size = 12),
    axis.text = ggplot2::element_text(size = 12),
    legend.position = "right",
    legend.text = ggplot2::element_text(face = "italic", size = 8),
    legend.title = ggplot2::element_text(size = 9),
    legend.key.height = grid::unit(0.4, "cm")
  ) +
  ggplot2::labs(x = NULL, y = "Cumulative sampled cells\n(log-scale)")

p2b <- ggplot2::ggplot(cs[year >= 1990], ggplot2::aes(year, cs_frac, group = species)) +
  ggplot2::geom_line(alpha = 0.4, colour = "#CC79A7") +
  ggplot2::geom_smooth(ggplot2::aes(group = 1), method = "loess",
                       se = TRUE, colour = "#7B3294", linewidth = 1) +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  theme_inv() + theme_bw() +  ggplot2::scale_x_continuous(
    limits = c(1990, 2025),
    breaks = xx_breaks
  ) +  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 12),
    axis.text = ggplot2::element_text(size = 12)
  ) + 
  ggplot2::labs(x = "Year", y = "Citizen-science share\nof annual records")

fig2 <- patchwork::wrap_plots(p2a, p2b, ncol = 1)
ggsave2("Fig2_sampling_effort.svg", fig2, 7, 7)

# =============================================================================
# Figure 3 -- naive vs all three effort-corrected range curves
# =============================================================================
# Same four series as Fig. 4/Table 3 (naive + the three corrections actually
# used for spread-rate estimation), not the separate relative-occupancy
# diagnostic (rel_occ is never converted into a spread-rate estimate; see
# 02_range_dynamics.R). This keeps Fig. 3 consistent with what the rest of
# the paper compares, rather than illustrating a fourth, unrelated metric.
long <- curves[, .(species, year, naive = aoo, rarefied = rar_aoo,
                   rarefied_regional = rarr_aoo, listlength = ll_aoo)] |>
  data.table::melt(id.vars = c("species", "year"),
                   variable.name = "method", value.name = "value")
# rescale each curve to 0-1 within species for comparability of shape
long[, value_s := value / max(value, na.rm = TRUE), by = .(species, method)]
long[, method := factor(method, levels = names(method_pal))]

fig3 <- ggplot2::ggplot(long[is.finite(value_s)],
                        ggplot2::aes(year, value_s, colour = method)) +
  ggplot2::geom_line(linewidth = 0.6, alpha = 0.85, na.rm = TRUE) +
  ggplot2::scale_colour_manual(values = method_pal, labels = method_lab, name = "Method") +
  ggplot2::facet_wrap(~species, ncol = 4, scales = "free_x") +
  theme_inv(9) +  theme_bw() +
  ggplot2::theme(strip.text = ggplot2::element_text(size = 12, face = "italic")) +
  ggplot2::labs(x = "Year", y = "Range metric (scaled 0-1 within species)") +  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 12),
    axis.text = ggplot2::element_text(size = 12)
  )  +  ggplot2::scale_x_continuous(
    limits = c(1950, 2025),
    breaks = x_breaks
  ) 

ggsave2("Fig3_range_curves.svg", fig3, 11, 8)

# =============================================================================
# Figure 4 -- spread velocity across all four correction methods (95% CI)
# =============================================================================
sr <- data.table::as.data.table(rd$spread)
sr[, species := factor(species, levels = sr[method == "naive"][order(rate_km_yr), species])]
sr[, method := factor(method, levels = names(method_pal))]

fig4 <- ggplot2::ggplot(sr, ggplot2::aes(rate_km_yr, species, colour = method)) +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = ci_lwr, xmax = ci_upr),
                          height = 0.25, position = ggplot2::position_dodge(0.6)) +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.6), size = 2) +
  ggplot2::scale_colour_manual(values = method_pal, labels = method_lab, name = "Method") +
  theme_inv() + theme_bw() +  
ggplot2::theme(
    axis.title = ggplot2::element_text(size = 12),
    axis.text = ggplot2::element_text(size = 12)
  ) +
  ggplot2::theme(axis.text.y = ggplot2::element_text(face = "italic")) +
  ggplot2::labs(x = expression("Radial spread rate (km yr"^-1*")"), y = NULL)
ggsave2("Fig4_spread_rates.svg", fig4, 9, 6.5)

# =============================================================================
# Figure 5 -- bias ratio vs citizen-science share, all three corrections
# =============================================================================
# The naive/rarefied ratio is the headline metric (Results 3.4), but the
# "does the bias scale with citizen-science reliance" test is repeated
# against all three corrections so the finding doesn't hinge on which
# correction is used as the denominator (cf. Fig. 6, which checks whether
# the corrections agree with each other on absolute spread rate).
options(scipen =9999)

bias <- data.table::as.data.table(rd$bias)
bias_long <- data.table::melt(
  bias[, .(species, group, cs_frac_recent, bias_ratio, bias_ratio_regional, bias_ratio_listlength)],
  id.vars = c("species", "group", "cs_frac_recent"),
  measure.vars = c("bias_ratio", "bias_ratio_regional", "bias_ratio_listlength"),
  variable.name = "method", value.name = "ratio")
bias_long[, method := factor(method,
  levels = c("bias_ratio", "bias_ratio_regional", "bias_ratio_listlength"),
  labels = c(method_lab[["rarefied"]], method_lab[["rarefied_regional"]], method_lab[["listlength"]]))]

# species with an undefined ratio (near-zero or NA denominator; safe_ratio()
# in 02_range_dynamics.R) have no meaningful y-position -- shown as rug ticks
# on the x-axis instead of being silently dropped, so it's visible *who* is
# missing from each panel and *where* they'd sit on the citizen-science axis.
bias_missing <- bias_long[!is.finite(ratio)]
bias_long <- bias_long[bias_long$stability !="Unstable (near-zero denominator)",]

fig5 <- ggplot2::ggplot(bias_long[is.finite(ratio)],
                        ggplot2::aes(cs_frac_recent, ratio, colour = group)) +
  ggplot2::geom_hline(yintercept = 1, linetype = 2, colour = "grey60") +
  ggplot2::geom_smooth(ggplot2::aes(group = 1), method = "lm", se = TRUE,
                       colour = "grey30", linewidth = 0.7) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_rug(data = bias_missing, ggplot2::aes(x = cs_frac_recent, colour = group),
                    sides = "b", inherit.aes = FALSE, length = grid::unit(0.06, "npc"),
                    linewidth = 0.8, show.legend = FALSE) +
  ggplot2::scale_colour_manual(values = group_pal, name = "Group") +
  ggplot2::scale_x_continuous(labels = scales::percent) +
  ggplot2::facet_wrap(~method, scales = "free_y") +
  theme_inv() + theme_bw() +
     ggplot2::theme(
    axis.title = ggplot2::element_text(size = 12),
    axis.text = ggplot2::element_text(size = 12)
  ) +
  ggplot2::labs(x = "Recent citizen-science share (2010-2024)",
                y = "Bias ratio (naive / corrected spread rate)",
                caption = "Tick marks on the x-axis = species with an undefined ratio (near-zero or NA denominator); shown at their citizen-science share only, since their ratio has no meaningful y-position.")
ggsave2("Fig5_bias_synthesis.svg", fig5, 13, 5)

# =============================================================================
# Figure 6 -- pairwise agreement among the three effort corrections
# (robustness check; naive is deliberately excluded -- it is the baseline
# being tested, not a candidate correction; see Figs. 4-5 for its bias)
# =============================================================================
corrected_methods <- c("rarefied", "rarefied_regional", "listlength")
method_pairs <- utils::combn(corrected_methods, 2, simplify = FALSE)

agree_long <- purrr::map_dfr(method_pairs, function(p) {
  data.table::data.table(
    species = bias$species, group = bias$group,
    rate_a = bias[[p[1]]], rate_b = bias[[p[2]]],
    pair = paste(method_lab[p[1]], "vs", method_lab[p[2]]))
})
agree_long <- data.table::as.data.table(agree_long)
agree_long[, pair := factor(pair, levels = unique(pair))]

fig6 <- ggplot2::ggplot(agree_long[is.finite(rate_a) & is.finite(rate_b)],
                        ggplot2::aes(rate_a, rate_b, colour = group)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::scale_colour_manual(values = group_pal, name = "Group") +
  ggplot2::facet_wrap(~pair, scales = "free", ncol = 3) +
  theme_inv() + theme_bw() +   ggplot2::theme(
    axis.title = ggplot2::element_text(size = 12),
    axis.text = ggplot2::element_text(size = 12)
  ) +
  ggplot2::labs(x = expression("Spread rate, method A (km yr"^-1*")"),
                y = expression("Spread rate, method B (km yr"^-1*")"))
ggsave2("Fig6_method_agreement.svg", fig6, 12, 4.5)

# =============================================================================
# Figure S2 -- leave-one-species-out sensitivity, rarefaction-based corrections
# =============================================================================
# The list-length version is omitted here: after excluding near-zero-
# denominator ratios (see bias_tbl / safe_ratio() in 02_range_dynamics.R),
# only 7 of 14 species retain a defined list-length-based bias ratio, and
# those 7 exclude every species with a citizen-science share above ~12% --
# too few and too range-restricted to jackknife informatively (Fig. 5 shows
# it as a scatter instead).
jk <- data.table::as.data.table(rd$jackknife_bias_cor)
jk <- jk[method %in% c("rarefied", "rarefied_regional")]
jk[, method := factor(method, levels = c("rarefied", "rarefied_regional"),
                      labels = method_lab[c("rarefied", "rarefied_regional")])]
# order species consistently across panels (by the headline/rarefied result)
order_levels <- jk[method == method_lab[["rarefied"]]][order(r), species_dropped]
jk[, species_dropped := factor(species_dropped, levels = order_levels)]

full_r_dt <- data.table::data.table(
  method = factor(method_lab[c("rarefied", "rarefied_regional")],
                  levels = levels(jk$method)),
  full_r = c(unname(rd$bias_test$estimate),
            unname(rd$bias_test_regional$estimate)))

figS2 <- ggplot2::ggplot(jk, ggplot2::aes(r, species_dropped)) +
  ggplot2::geom_vline(data = full_r_dt, ggplot2::aes(xintercept = full_r),
                      linetype = 2, colour = "#D55E00") +
  ggplot2::geom_vline(xintercept = 0, linetype = 3, colour = "grey60") +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = ci_lwr, xmax = ci_upr), height = 0.2) +
  ggplot2::geom_point(size = 2, colour = "#0072B2") +
  ggplot2::facet_wrap(~method, scales = "free_x") +
  theme_inv() +
  ggplot2::theme(axis.text.y = ggplot2::element_text(face = "italic")) +
  ggplot2::labs(x = expression("Pearson "*italic(r)*" (bias ratio ~ citizen-science share)"),
                y = "Species excluded",
                title = "Leave-one-species-out sensitivity, rarefaction-based corrections",
                caption = "Dashed line = full-sample estimate for that correction; dotted line = zero.")
ggsave2("FigS2_jackknife_bias_cor.svg", figS2, 9, 6)

# =============================================================================
# Tables
# =============================================================================
n_summary <- purrr::imap(occ$focal, ~ data.table::data.table(
  species = .y, n_records = nrow(.x),
  year_first = min(.x$year), year_last = max(.x$year))) |>
  data.table::rbindlist()

doi_tbl <- tryCatch(data.table::fread(file.path(cfg$dir_tab, "gbif_download_dois.csv")),
                    error = function(e) data.table::data.table())

table1 <- sp_t |>
  dplyr::select(species, group, native_region) |>
  dplyr::left_join(n_summary, by = "species")
readr::write_csv(table1, file.path(cfg$dir_tab, "Table1_species.csv"))

table2 <- sr |>
  dplyr::select(species, group, method, lag_end, sat_start) |>
  tidyr::pivot_wider(names_from = method, values_from = c(lag_end, sat_start))
readr::write_csv(table2, file.path(cfg$dir_tab, "Table2_phases.csv"))

table3 <- sr |>
  dplyr::transmute(species, group, method,
                   rate = sprintf("%.1f (%.1f-%.1f)", rate_km_yr, ci_lwr, ci_upr)) |>
  tidyr::pivot_wider(names_from = method, values_from = rate)
readr::write_csv(table3, file.path(cfg$dir_tab, "Table3_spread.csv"))

# =============================================================================
# results_values.json  --  numbers injected into the manuscript template
# =============================================================================
med_naive <- stats::median(sr[method == "naive", rate_km_yr], na.rm = TRUE)
med_rar   <- stats::median(sr[method == "rarefied", rate_km_yr], na.rm = TRUE)
med_rarr  <- stats::median(sr[method == "rarefied_regional", rate_km_yr], na.rm = TRUE)
med_ll    <- stats::median(sr[method == "listlength", rate_km_yr], na.rm = TRUE)

# agreement among the three *corrected* methods (excludes naive)
corrected_methods <- c("rarefied", "rarefied_regional", "listlength")
ma_corrected <- dplyr::filter(rd$method_agreement,
                              method_a %in% corrected_methods, method_b %in% corrected_methods)

# leave-one-species-out sensitivity of the bias~citizen-science correlation
# (jkc = headline/continental-rarefaction jackknife, used throughout Results
# 3.4; jkc_regional/jkc_listlength are the added checks for the other two
# corrections, summarized only by their r range, not per-species detail)
jkc            <- dplyr::filter(rd$jackknife_bias_cor, method == "rarefied")
jkc_regional   <- dplyr::filter(rd$jackknife_bias_cor, method == "rarefied_regional")
jkc_listlength <- dplyr::filter(rd$jackknife_bias_cor, method == "listlength")
jk_r_drop <- function(sp) round(jkc$r[jkc$species_dropped == sp], 2)

res_vals <- list(
  n_species          = nrow(sp_t),
  n_records_total    = sum(n_summary$n_records),
  n_records_median   = as.integer(stats::median(n_summary$n_records)),
  year_range         = paste(cfg$year_min, cfg$year_max, sep = "-"),
  grid_km            = cfg$cell_km,
  accessible_buffer_km    = cfg$accessible_buffer_km,
  median_spread_naive     = round(med_naive, 1),
  median_spread_rarefied  = round(med_rar, 1),
  median_spread_regional  = round(med_rarr, 1),
  median_spread_listlength= round(med_ll, 1),
  median_bias_ratio  = round(stats::median(bias$bias_ratio[is.finite(bias$bias_ratio)]), 2),
  bias_cor_r         = round(unname(rd$bias_test$estimate), 2),
  bias_cor_ci        = round(as.numeric(rd$bias_test$conf.int), 2),
  bias_cor_p         = signif(rd$bias_test$p.value, 2),
  bias_cor_r_regional    = round(unname(rd$bias_test_regional$estimate), 2),
  bias_cor_p_regional    = signif(rd$bias_test_regional$p.value, 2),
  bias_cor_r_listlength  = round(unname(rd$bias_test_listlength$estimate), 2),
  bias_cor_p_listlength  = signif(rd$bias_test_listlength$p.value, 2),
  kendall_w              = round(rd$kendall_w, 2),
  agreement_rho_min      = round(min(ma_corrected$rho, na.rm = TRUE), 2),
  agreement_rho_max      = round(max(ma_corrected$rho, na.rm = TRUE), 2),
  jackknife_r_min             = round(min(jkc$r), 2),
  jackknife_r_max             = round(max(jkc$r), 2),
  jackknife_p_min             = signif(min(jkc$p_value), 2),
  jackknife_r_drop_dbugensis  = jk_r_drop("Dreissena bugensis"),
  jackknife_r_drop_amelas     = jk_r_drop("Ameiurus melas"),
  jackknife_r_min_regional    = round(min(jkc_regional$r), 2),
  jackknife_r_max_regional    = round(max(jkc_regional$r), 2),
  jackknife_r_min_listlength  = round(min(jkc_listlength$r), 2),
  jackknife_r_max_listlength  = round(max(jkc_listlength$r), 2),
  fastest_species    = sr[method == "rarefied"][which.max(rate_km_yr), as.character(species)],
  generated          = as.character(Sys.Date())
)
jsonlite::write_json(res_vals, file.path(cfg$dir_out, "results_values.json"),
                     auto_unbox = TRUE, pretty = TRUE)

message("Figures + tables written to ", cfg$dir_fig, " and ", cfg$dir_tab)
print(res_vals)
