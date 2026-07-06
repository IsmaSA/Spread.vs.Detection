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
# Figure 1 -- occurrence density across Europe
# =============================================================================
map_bbox <- sf::st_bbox(c(xmin = 2.5e6, xmax = 6e6, ymin = 1.5e6, ymax = 5.4e6),
                        crs = cfg$crs_analysis)

europe <- rnaturalearth::ne_countries(scale = "medium", continent = "Europe",
                                       returnclass = "sf") |>
  sf::st_transform(cfg$crs_analysis)

pts <- purrr::imap(occ$focal, ~ .x[, .(lon, lat)][, species := .y]) |>
  data.table::rbindlist() |>
  sf::st_as_sf(coords = c("lon", "lat"), crs = cfg$crs_geo) |>
  sf::st_transform(cfg$crs_analysis)
pts <- dplyr::left_join(pts, dplyr::select(sp_t, species, group), by = "species")
# randomize draw order across groups so no single group is systematically
# overplotted on top of the others in the pooled panel (an easy-to-miss bias
# when point layers are drawn in a fixed, data-order sequence)
pts_pooled <- pts[sample.int(nrow(pts)), ]

# ---- panel (a): all species pooled, coloured by taxonomic group -----------
fig1a <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = europe, fill = "grey96", colour = "grey80", linewidth = 0.2) +
  ggplot2::geom_sf(data = pts_pooled, ggplot2::aes(colour = group), size = 0.15, alpha = 0.3) +
  ggplot2::scale_colour_manual(values = group_pal, name = "Group") +
  ggplot2::coord_sf(crs = sf::st_crs(cfg$crs_analysis),
                    xlim = map_bbox[c("xmin", "xmax")], ylim = map_bbox[c("ymin", "ymax")]) +
  ggspatial::annotation_scale(location = "bl", height = grid::unit(0.15, "cm")) +
  theme_inv(10) +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1))) +
  ggplot2::labs(subtitle = "(a) All 14 species pooled, by taxonomic group")

# ---- panel (b): faceted by species -----------------------------------------
fig1b <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = europe, fill = "grey96", colour = "grey80", linewidth = 0.2) +
  ggplot2::geom_sf(data = pts, ggplot2::aes(colour = group), size = 0.15, alpha = 0.35) +
  ggplot2::scale_colour_manual(values = group_pal, name = "Group") +
  ggplot2::coord_sf(crs = sf::st_crs(cfg$crs_analysis),
                    xlim = map_bbox[c("xmin", "xmax")], ylim = map_bbox[c("ymin", "ymax")]) +
  ggplot2::facet_wrap(~species, ncol = 4) +
  theme_inv(9) +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2, alpha = 1))) +
  ggplot2::labs(subtitle = "(b) Faceted by species")

fig1 <- patchwork::wrap_plots(fig1a, fig1b, ncol = 1, heights = c(0.3, 0.7)) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title = "Cleaned GBIF occurrences of study invaders (1950-2024)",
    caption = "ETRS89-LAEA (EPSG:3035). Points = individual records.")

ggsave2("Fig1_occurrence_maps.svg", fig1, 11, 13)

# =============================================================================
# Figure 2 -- sampling effort through time + citizen-science share
# =============================================================================
effort_ts <- curves[, .(sampled = max(sampled, na.rm = TRUE)), by = .(species, year)]
cs <- data.table::fread(file.path(cfg$dir_tab, "citizen_science_share.csv"))

p2a <- ggplot2::ggplot(effort_ts, ggplot2::aes(year, sampled, group = species)) +
  ggplot2::geom_line(alpha = 0.5, colour = "#0072B2") +
  ggplot2::scale_y_log10() +
  theme_inv() +
  ggplot2::labs(x = NULL, y = "Cumulative sampled cells\n(target-group footprint, log)",
                subtitle = "(a) Sampling-effort footprint accumulates over time")

p2b <- ggplot2::ggplot(cs[year >= 1990], ggplot2::aes(year, cs_frac, group = species)) +
  ggplot2::geom_line(alpha = 0.4, colour = "#CC79A7") +
  ggplot2::geom_smooth(ggplot2::aes(group = 1), method = "loess",
                       se = TRUE, colour = "#7B3294", linewidth = 1) +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  theme_inv() +
  ggplot2::labs(x = "Year", y = "Citizen-science share\nof annual records",
                subtitle = "(b) Rising share of citizen-science (iNaturalist) records")

fig2 <- patchwork::wrap_plots(p2a, p2b, ncol = 1)
ggsave2("Fig2_sampling_effort.png", fig2, 7, 7)

# =============================================================================
# Figure 3 -- naive vs effort-corrected range curves
# =============================================================================
long <- curves[, .(species, year, Naive = aoo, `Relative occupancy` = rel_occ)] |>
  data.table::melt(id.vars = c("species", "year"),
                   variable.name = "metric", value.name = "value")
# rescale each metric to 0-1 within species for comparability of shape
long[, value_s := value / max(value, na.rm = TRUE), by = .(species, metric)]

fig3 <- ggplot2::ggplot(long[is.finite(value_s)],
                        ggplot2::aes(year, value_s, colour = metric)) +
  ggplot2::geom_line(linewidth = 0.7, na.rm = TRUE) +
  ggplot2::scale_colour_manual(values = c("Naive" = "#D55E00",
                                          "Relative occupancy" = "#0072B2"),
                               name = NULL) +
  ggplot2::facet_wrap(~species, ncol = 4, scales = "free_x") +
  theme_inv(9) +
  ggplot2::labs(x = "Year", y = "Range metric (scaled 0-1 within species)",
                title = "Naive range accumulation vs effort-corrected occupancy",
                caption = paste("Naive = cumulative first-detection cells;",
                                "Relative occupancy = occupied / sampled cells."))
ggsave2("Fig3_range_curves.png", fig3, 11, 8)

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
  theme_inv() +
  ggplot2::theme(axis.text.y = ggplot2::element_text(face = "italic")) +
  ggplot2::labs(x = expression("Radial spread rate (km yr"^-1*")"), y = NULL,
                title = "Effort correction lowers estimated spread velocity",
                subtitle = "Naive vs. three correction methods (rarefaction, accessible-area rarefaction, list-length model)")
ggsave2("Fig4_spread_rates.png", fig4, 9, 6.5)

# =============================================================================
# Figure 5 -- bias ratio vs citizen-science share (synthesis)
# =============================================================================
bias <- data.table::as.data.table(rd$bias)
fig5 <- ggplot2::ggplot(bias[is.finite(bias_ratio)],
                        ggplot2::aes(cs_frac_recent, bias_ratio, colour = group)) +
  ggplot2::geom_hline(yintercept = 1, linetype = 2, colour = "grey60") +
  ggplot2::geom_smooth(ggplot2::aes(group = 1), method = "lm", se = TRUE,
                       colour = "grey30", linewidth = 0.7) +
  ggplot2::geom_point(size = 3) +
  ggplot2::scale_colour_manual(values = group_pal, name = "Group") +
  ggplot2::scale_x_continuous(labels = scales::percent) +
  theme_inv() +
  ggplot2::labs(x = "Recent citizen-science share (2010-2024)",
                y = "Bias ratio (naive / corrected spread rate)",
                title = "Naive spread inflation scales with citizen-science share")
ggsave2("Fig5_bias_synthesis.png", fig5, 7, 6)

# =============================================================================
# Figure 6 -- agreement among correction methods (robustness check)
# =============================================================================
agree_long <- data.table::melt(
  bias[, .(species, group, rarefied, rarefied_regional, listlength, naive)],
  id.vars = c("species", "group", "rarefied"),
  measure.vars = c("rarefied_regional", "listlength", "naive"),
  variable.name = "comparator", value.name = "rate_comparator")
agree_long[, comparator := factor(comparator,
  levels = c("rarefied_regional", "listlength", "naive"),
  labels = c("Rarefied (accessible area)", "List-length model", "Naive"))]

fig6 <- ggplot2::ggplot(agree_long[is.finite(rarefied) & is.finite(rate_comparator)],
                        ggplot2::aes(rarefied, rate_comparator, colour = group)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::scale_colour_manual(values = group_pal, name = "Group") +
  ggplot2::facet_wrap(~comparator, scales = "free") +
  theme_inv() +
  ggplot2::labs(x = expression("Spread rate, rarefied (continental), km yr"^-1),
                y = expression("Spread rate, alternative method, km yr"^-1),
                title = "Agreement among correction methods",
                caption = "Dashed line = 1:1 agreement with the primary (continental rarefaction) estimate.")
ggsave2("Fig6_method_agreement.png", fig6, 10, 4)

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
  kendall_w              = round(rd$kendall_w, 2),
  agreement_rho_min      = round(min(ma_corrected$rho, na.rm = TRUE), 2),
  agreement_rho_max      = round(max(ma_corrected$rho, na.rm = TRUE), 2),
  fastest_species    = sr[method == "rarefied"][which.max(rate_km_yr), as.character(species)],
  generated          = as.character(Sys.Date())
)
jsonlite::write_json(res_vals, file.path(cfg$dir_out, "results_values.json"),
                     auto_unbox = TRUE, pretty = TRUE)

message("Figures + tables written to ", cfg$dir_fig, " and ", cfg$dir_tab)
print(res_vals)
