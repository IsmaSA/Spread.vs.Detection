# =============================================================================
# 03_synthesis_figures_tables.R  --  publication-ready figures + tables
#
# Figures (300 dpi, colourblind-safe viridis/Okabe-Ito, EPSG:3035 maps):
#   Fig 1  occurrence density of study species across Europe
#   Fig 2  sampling effort through time (the confound) + citizen-science share
#   Fig 3  naive vs effort-corrected range curves per species
#   Fig 4  spread velocity, naive vs rarefied, with 95% CI
#   Fig 5  cross-species synthesis: bias ratio vs citizen-science share
# Tables:
#   Table 1  species, group, native region, n records, span, DOI
#   Table 2  invasion-phase breakpoints (naive vs corrected)
#   Table 3  spread rates (km/yr) with CIs
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
europe <- rnaturalearth::ne_countries(scale = "medium", continent = "Europe",
                                       returnclass = "sf") |>
  sf::st_transform(cfg$crs_analysis)

pts <- purrr::imap(occ$focal, ~ .x[, .(lon, lat)][, species := .y]) |>
  data.table::rbindlist() |>
  sf::st_as_sf(coords = c("lon", "lat"), crs = cfg$crs_geo) |>
  sf::st_transform(cfg$crs_analysis)
pts <- dplyr::left_join(pts, dplyr::select(sp_t, species, group), by = "species")

fig1 <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = europe, fill = "grey96", colour = "grey80", linewidth = 0.2) +
  ggplot2::geom_sf(data = pts, ggplot2::aes(colour = group), size = 0.25, alpha = 0.35) +
  ggplot2::scale_colour_manual(values = group_pal, name = "Group") +
  ggplot2::coord_sf(crs = sf::st_crs(cfg$crs_analysis),
                    xlim = c(2.5e6, 6e6), ylim = c(1.5e6, 5.4e6)) +
  ggspatial::annotation_scale(location = "bl", height = grid::unit(0.15, "cm")) +
  ggplot2::facet_wrap(~species, ncol = 4) +
  theme_inv(9) +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2, alpha = 1))) +
  ggplot2::labs(title = "Cleaned GBIF occurrences of study invaders (1950-2024)",
                caption = "ETRS89-LAEA (EPSG:3035). Points = individual records.")
ggsave2("Fig1_occurrence_maps.png", fig1, 11, 9)

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
# Figure 4 -- spread velocity naive vs rarefied (95% CI)
# =============================================================================
sr <- data.table::as.data.table(rd$spread)
sr[, species := factor(species, levels = sr[method == "naive"][order(rate_km_yr), species])]

fig4 <- ggplot2::ggplot(sr, ggplot2::aes(rate_km_yr, species, colour = method)) +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = ci_lwr, xmax = ci_upr),
                          height = 0.25, position = ggplot2::position_dodge(0.5)) +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5), size = 2) +
  ggplot2::scale_colour_manual(values = c(naive = "#D55E00", rarefied = "#0072B2"),
                               name = "Estimate") +
  theme_inv() +
  ggplot2::theme(axis.text.y = ggplot2::element_text(face = "italic")) +
  ggplot2::labs(x = expression("Radial spread rate (km yr"^-1*")"), y = NULL,
                title = "Effort correction lowers estimated spread velocity")
ggsave2("Fig4_spread_rates.png", fig4, 8, 6)

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
res_vals <- list(
  n_species          = nrow(sp_t),
  n_records_total    = sum(n_summary$n_records),
  n_records_median   = as.integer(stats::median(n_summary$n_records)),
  year_range         = paste(cfg$year_min, cfg$year_max, sep = "-"),
  grid_km            = cfg$cell_km,
  median_spread_naive    = round(med_naive, 1),
  median_spread_rarefied = round(med_rar, 1),
  median_bias_ratio  = round(stats::median(bias$bias_ratio[is.finite(bias$bias_ratio)]), 2),
  bias_cor_r         = round(unname(rd$bias_test$estimate), 2),
  bias_cor_ci        = round(as.numeric(rd$bias_test$conf.int), 2),
  bias_cor_p         = signif(rd$bias_test$p.value, 2),
  fastest_species    = sr[method == "rarefied"][which.max(rate_km_yr), as.character(species)],
  generated          = as.character(Sys.Date())
)
jsonlite::write_json(res_vals, file.path(cfg$dir_out, "results_values.json"),
                     auto_unbox = TRUE, pretty = TRUE)

message("Figures + tables written to ", cfg$dir_fig, " and ", cfg$dir_tab)
print(res_vals)
