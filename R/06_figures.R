# Figures 1-6 and S1-S5.
source("R/functions.R")
library(ggplot2)
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)
tb <- function(f) fread(file.path("output/tables", f))
save_fig <- function(name, p, w, h) {
  ggsave(file.path("output/figures", paste0(name, ".png")), p, width = w, height = h, dpi = 300, device = ragg::agg_png, bg = "white")
  ggsave(file.path("output/figures", paste0(name, ".svg")), p, width = w, height = h, device = svglite::svglite)
}
grp_pal <- c(Crustacean = "#B2182B", Mollusc = "#6A3D9A", Fish = "#4D4D4D", Plant = "#E6AB02")
method_pal <- c(naive = "#D55E00", rarefied = "#0072B2", accessible = "#009E73", listlength = "#CC79A7")
mlab <- c(naive = "Naïve", rarefied = "Rarefied (continental)", accessible = "Rarefied (accessible area)", listlength = "List-length model")
dag <- function(s) ifelse(sp$intro_year[match(s, sp$species)] < 1950, paste0(s, " †"), s)
abbr <- function(s) sub("^([A-Z])[a-z]+ ", "\\1. ", s)

ne <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |> st_transform(3035)
box <- native_box()

# Figure 1
eu <- ne[ne$continent == "Europe", ]
tk <- ne[ne$name == "Turkey" | ne$admin == "Turkey", ]
cells <- rbindlist(lapply(sp$species, function(s) {
  g <- grid(readRDS(taxon_file(s)), 25)
  unique(g[, .(gx, gy)])[, species := s]
}))
cells[, `:=`(cx = (gx + 0.5) * 25000, cy = (gy + 0.5) * 25000)]
cells <- merge(cells, sp[, .(species, group)], by = "species")
cells[, species := factor(species, levels = sort(unique(species)))]
boxes <- st_sf(species = factor(pc, levels = levels(cells$species)), geometry = rep(box, length(pc)))
p <- ggplot() +
  geom_sf(data = eu, fill = "grey88", colour = "grey60", linewidth = 0.1) +
  geom_sf(data = tk, fill = "grey96", colour = "grey75", linewidth = 0.1) +
  geom_tile(data = cells, aes(cx, cy, fill = group), width = 25000, height = 25000) +
  geom_sf(data = boxes, fill = NA, colour = "black", linewidth = 0.3) +
  scale_fill_manual(values = grp_pal, name = "Taxonomic group") +
  coord_sf(crs = st_crs(3035), xlim = c(2.5e6, 6.2e6), ylim = c(1.4e6, 5.4e6), datum = NA) +
  facet_wrap(~species, ncol = 4) +
  theme_void(base_size = 8) +
  theme(strip.text = element_text(face = "italic", size = 8), legend.position = "bottom",
        panel.border = element_rect(colour = "grey60", fill = NA, linewidth = 0.3))
save_fig("fig1_maps", p, 6.8, 7.6)

# Figure 2
eff <- merge(tb("sampled_cells.csv"), unique(sp[, .(target, group)]), by = "target")
lt <- c(Bivalvia = "solid", Cambaridae = "solid", Astacidae = "dashed", Amphipoda = "dotted", Perciformes = "solid",
        Cypriniformes = "dashed", Cyprinodontiformes = "dotted", Siluriformes = "dotdash", Hydrocharitaceae = "solid", Onagraceae = "dashed")
txt <- theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_text(size = 13))
pa <- ggplot(eff[cells > 0], aes(year, cells, colour = group, linetype = target)) +
  geom_line(linewidth = 0.6) +
  scale_y_log10(labels = scales::label_comma()) +
  scale_colour_manual(values = grp_pal, name = "Taxonomic group") +
  scale_linetype_manual(values = lt, name = "Target group") +
  scale_x_continuous(limits = c(1950, 2025), breaks = seq(1950, 2020, 10)) +
  labs(x = NULL, y = "Cumulative cells sampled", tag = "a") +
  theme_bw(base_size = 11) + txt +
  theme(legend.key.width = grid::unit(1, "cm"), legend.key.height = grid::unit(0.38, "cm"), legend.spacing.y = grid::unit(0.1, "cm"),
        legend.text = element_text(size = 9), legend.title = element_text(size = 10), legend.margin = margin(0, 0, 0, 0))
cs <- merge(tb("inat_share_by_year.csv")[year >= 1990], sp[, .(species, group)], by = "species")
hl <- data.table(species = c("Ameiurus melas", "Gambusia holbrooki", "Dreissena polymorpha", "Elodea nuttallii"),
                 lx = 2013, ly = c(0.88, 0.71, 0.55, 0.14))
hl <- merge(hl, cs[year == 2024, .(species, group, ey = inat)], by = "species")
pb <- ggplot(cs, aes(year, inat, group = species, colour = group)) +
  geom_line(alpha = 0.8, linewidth = 0.45) +
  geom_segment(data = hl, aes(x = lx + 0.5, y = ly + 0.015, xend = 2023.6, yend = ey, colour = group), inherit.aes = FALSE,
               arrow = grid::arrow(length = grid::unit(0.18, "cm"), type = "closed"), linewidth = 0.4, show.legend = FALSE) +
  geom_text(data = hl, aes(x = lx, y = ly, label = species, colour = group), inherit.aes = FALSE, hjust = 1, vjust = 0,
            fontface = "italic", size = 3.6, show.legend = FALSE) +
  scale_colour_manual(values = grp_pal) +
  scale_y_continuous(labels = scales::label_percent()) +
  scale_x_continuous(limits = c(1990, 2025), breaks = seq(1990, 2020, 10)) +
  labs(x = "Year", y = "iNaturalist share of\nannual occurrences", tag = "b") +
  theme_bw(base_size = 11) + txt + theme(legend.position = "none")
save_fig("fig2_effort_inat", patchwork::wrap_plots(pa, pb, ncol = 1), 7.5, 6.8)

# Figure 3
cv <- tb("range_curves_25km.csv")[method %in% meth]
cv[, value := aoo / max(aoo), by = .(species, method)]
cv[, method := factor(method, levels = meth)]
cv[, lab := factor(dag(species), levels = dag(sort(unique(species))))]
p <- ggplot(cv[is.finite(value)], aes(year, value, colour = method)) +
  geom_line(linewidth = 0.45) +
  scale_colour_manual(values = method_pal, labels = mlab, name = "Method") +
  scale_x_continuous(limits = c(1950, 2025), breaks = c(1950, 1980, 2010)) +
  facet_wrap(~lab, ncol = 4) +
  labs(x = "Year", y = "Occupied cells") +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(face = "italic", size = 9), axis.text = element_text(size = 9), axis.title = element_text(size = 11),
        legend.position = "bottom", legend.text = element_text(size = 10), legend.title = element_text(size = 10))
save_fig("fig3_range_curves", p, 7, 6.8)

# Figure 4
sr <- tb("spread_rates.csv")[km == 25 & method %in% meth]
clip <- 30
sr[, `:=`(x = pmin(rate, clip), lo = pmax(pmin(ci_lwr, clip), -3), hi = pmin(ci_upr, clip))]
sr[, species := factor(species, levels = sr[method == "naive"][order(rate), species])]
sr[, method := factor(method, levels = meth)]
p <- ggplot(sr, aes(x, species, colour = method)) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(0.7), size = 0.28, linewidth = 0.35) +
  scale_colour_manual(values = method_pal, labels = mlab, name = NULL, guide = guide_legend(nrow = 2, byrow = TRUE)) +
  coord_cartesian(xlim = c(-3, clip)) +
  labs(x = expression("Radial spread rate (km yr"^-1*")"), y = NULL) +
  theme_bw(base_size = 10) +
  theme(axis.text.y = element_text(face = "italic", size = 10), axis.text.x = element_text(size = 10), axis.title = element_text(size = 11),
        legend.position = "bottom", legend.text = element_text(size = 10), legend.key.width = grid::unit(0.8, "cm"),
        legend.spacing.x = grid::unit(0.15, "cm"), legend.margin = margin(0, 0, 0, 0))
save_fig("fig4_spread_rates", p, 7.3, 5.8)

# Figure 5
ss <- tb("simulation_runs.csv")[is.finite(ratio) & method %in% meth]
ss[, sc := factor(scenario, levels = c("effort", "boost"), labels = c("Effort scenario", "Boost scenario"))]
ss[, ml := factor(method, levels = rev(meth), labels = rev(mlab))]
med <- ss[, .(med = median(ratio)), by = .(sc, method, ml)]
brk <- c(0.05, 0.1, 0.25, 0.5, 1, 2, 4, 10, 20)
p <- ggplot(ss, aes(ratio, ml, colour = method)) +
  geom_rect(data = data.frame(sc = factor(levels(ss$sc), levels = levels(ss$sc))), aes(xmin = 0.75, xmax = 1.25, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "grey85", alpha = 0.6) +
  geom_vline(xintercept = 1, colour = "grey40", linewidth = 0.4) +
  geom_point(position = position_jitter(height = 0.2, width = 0, seed = 1), size = 1.5, alpha = 0.45) +
  geom_segment(data = med, aes(x = 1, xend = med, yend = ml), linewidth = 1.1) +
  geom_point(data = med, aes(x = med), size = 4.2) +
  geom_text(data = med, aes(x = med, label = sprintf("%.2f", med)), vjust = -1.5, size = 3.2, colour = "grey20", show.legend = FALSE) +
  scale_x_log10(breaks = brk, labels = as.character(brk)) +
  coord_cartesian(xlim = c(0.03, 25)) +
  scale_colour_manual(values = method_pal, guide = "none") +
  facet_wrap(~sc, ncol = 1) +
  labs(x = "Estimates ⁄ true spread rate", y = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 9.5),
        axis.title = element_text(size = 11), strip.text = element_text(size = 10))
save_fig("fig5_simulation", p, 6.4, 5.8)

# Figure 6
b <- tb("reduction_and_shares.csv")
cr <- tb("correlations.csv")
corr <- meth[-1]
bl <- melt(b[, c("species", "group", "inat", paste0("red_", corr)), with = FALSE], id.vars = c("species", "group", "inat"),
           variable.name = "method", value.name = "red")
bl[, method := factor(sub("red_", "", method), levels = corr)]
ann <- cr[predictor == "inat" & response %in% paste0("red_", corr),
          .(method = factor(sub("red_", "", response), levels = corr), txt = sprintf("r = %.2f (%.2f to %.2f), P = %.2f", r, lwr, upr, p))]
p <- ggplot(bl, aes(inat, red, colour = group)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, se = TRUE, colour = "grey30", fill = "grey85", linewidth = 0.5) +
  geom_point(size = 1.8) +
  geom_text(data = ann, aes(x = Inf, y = Inf, label = txt), inherit.aes = FALSE, hjust = 1.04, vjust = 1.5, size = 2.6) +
  scale_colour_manual(values = grp_pal, name = "Taxonomic group") +
  scale_x_continuous(labels = scales::label_percent()) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  facet_wrap(~method, scales = "free_y", labeller = labeller(method = mlab)) +
  labs(x = "Mean annual iNaturalist share of occurrences", y = "Relative reduction in spread rate") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", strip.text = element_text(size = 9), legend.text = element_text(size = 9), legend.title = element_text(size = 9))
save_fig("fig6_reduction_inat", p, 6.8, 3.4)

# Figure S1
nb <- rbindlist(lapply(pc, function(s) {
  d <- readRDS(taxon_file(s))
  g <- unique(grid(d, 25)[, .(gx, gy)])
  g[, `:=`(cx = (gx + 0.5) * 25000, cy = (gy + 0.5) * 25000, species = s)]
  inb <- lengths(st_intersects(st_as_sf(g, coords = c("cx", "cy"), crs = 3035, remove = FALSE), box)) > 0
  g[, status := ifelse(inb, "Removed (inside box)", "Retained")]
  g[, lab := sprintf("%s\n(%d of %s records removed)", s, sum(d$box), format(nrow(d), big.mark = ","))]
}))
nb[, lab := factor(lab, levels = unique(lab))]
p <- ggplot() +
  geom_sf(data = ne, fill = "grey90", colour = "grey65", linewidth = 0.1) +
  geom_tile(data = nb, aes(cx, cy, fill = status), width = 25000, height = 25000) +
  geom_sf(data = box, fill = NA, colour = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("Removed (inside box)" = "#B2182B", Retained = "#2166AC"), name = NULL) +
  coord_sf(crs = st_crs(3035), xlim = c(2.6e6, 7.4e6), ylim = c(1.4e6, 5.6e6), datum = NA) +
  facet_wrap(~lab, ncol = 2) +
  theme_void(base_size = 8) +
  theme(strip.text = element_text(face = "italic", size = 8), legend.position = "bottom",
        panel.border = element_rect(colour = "grey60", fill = NA, linewidth = 0.3))
save_fig("figS1_native_box", p, 6.5, 6.2)

# Figure S2
fits <- rbindlist(lapply(meth, function(m) rbindlist(lapply(sp$species, function(s) {
  d <- cv[method == m & species == s & aoo > 0, .(year, aoo)]
  bp <- breaks(d)
  if (anyNA(bp)) return(NULL)
  X <- cbind(1, d$year, pmax(0, d$year - bp[1]), pmax(0, d$year - bp[2]))
  y <- log1p(d$aoo)
  data.table(species = s, method = m, year = d$year, y = y, fit = drop(X %*% lm.fit(X, y)$coefficients), bp1 = bp[1], bp2 = bp[2])
}))))
fits[, method := factor(method, levels = meth, labels = mlab)]
lv <- sort(unique(fits$species))
fits[, lab := factor(paste0(abbr(species), ifelse(sp$intro_year[match(species, sp$species)] < 1950, " †", "")),
                     levels = paste0(abbr(lv), ifelse(sp$intro_year[match(lv, sp$species)] < 1950, " †", "")))]
k <- c("Observed", "Fitted two-breakpoint regression", "End of lag phase", "Onset of saturation")
p <- ggplot(fits, aes(year)) +
  geom_point(aes(y = y, colour = k[1]), size = 0.5) +
  geom_line(aes(y = fit, colour = k[2]), linewidth = 0.6) +
  geom_vline(aes(xintercept = bp1, colour = k[3]), linetype = 2, linewidth = 0.35) +
  geom_vline(aes(xintercept = bp2, colour = k[4]), linetype = 3, linewidth = 0.45) +
  scale_colour_manual(name = NULL, breaks = k, values = setNames(c("grey35", "#D55E00", "black", "black"), k),
                      guide = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(shape = c(16, NA, NA, NA), linetype = c(0, 1, 2, 3),
                                                                                        linewidth = c(0, 0.7, 0.5, 0.6), size = c(1.4, NA, NA, NA)))) +
  facet_grid(lab ~ method, scales = "free", switch = "y") +
  scale_x_continuous(breaks = scales::breaks_pretty(3)) +
  scale_y_continuous(n.breaks = 3) +
  labs(x = "Year", y = "log(1 + occupied cells)") +
  theme_bw(base_size = 9) +
  theme(strip.text.y.left = element_text(face = "italic", angle = 0, size = 9, hjust = 1), strip.text.x = element_text(size = 8, face = "bold"),
        strip.background = element_blank(), strip.placement = "outside", axis.text = element_text(size = 7), axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 7, margin = margin(t = 1)), panel.spacing = grid::unit(0.18, "lines"), legend.position = "bottom",
        legend.text = element_text(size = 9), legend.key.width = grid::unit(1.1, "cm"), plot.margin = margin(3, 5, 3, 3))
save_fig("figS2_breakpoints", p, 7.6, 9.6)

# Figure S3
cv[, t := year - sp$intro_year[match(species, sp$species)]]
cv[, start := 1950 - sp$intro_year[match(species, sp$species)]]
l0 <- "Start of study period (1950)"
p <- ggplot(cv[is.finite(value)], aes(t, value, colour = method)) +
  geom_line(linewidth = 0.5) +
  geom_vline(aes(xintercept = start, colour = l0), linetype = 3, linewidth = 0.4) +
  scale_colour_manual(name = NULL, breaks = c(meth, l0), labels = c(mlab, l0), values = c(method_pal, setNames("grey40", l0)),
                      guide = guide_legend(ncol = 1, override.aes = list(linetype = c(1, 1, 1, 1, 3), linewidth = 0.7))) +
  facet_wrap(~lab, ncol = 4, scales = "free_x") +
  labs(x = "Years since first introduction", y = "Scaled curve") +
  theme_bw(base_size = 9) +
  theme(strip.text = element_text(face = "italic", size = 8), axis.text = element_text(size = 8), axis.title = element_text(size = 10),
        legend.position = "inside", legend.position.inside = c(0.765, 0.13), legend.justification.inside = c(0.5, 0.5),
        legend.text = element_text(size = 8.5), legend.key.width = grid::unit(0.9, "cm"), legend.margin = margin(4, 6, 4, 4),
        legend.background = element_rect(fill = "white", colour = "grey60", linewidth = 0.3))
save_fig("figS3_time_since_introduction", p, 7.4, 6.6)

# Figure S4
pr <- list(c("rarefied", "accessible"), c("rarefied", "listlength"), c("accessible", "listlength"))
pt <- paste0(LETTERS[1:3], ": ", c("Continental vs accessible-area rarefaction", "Continental rarefaction vs list-length model",
                                   "Accessible-area rarefaction vs list-length model"))
ag <- rbindlist(lapply(1:3, function(i) data.table(group = b$group, a = b[[pr[[i]][1]]], b = b[[pr[[i]][2]]], panel = factor(pt[i], levels = pt))))
pl <- scales::pseudo_log_trans(base = 10)
p <- ggplot(ag, aes(a, b, colour = group)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
  geom_point(size = 1.8) +
  scale_colour_manual(values = grp_pal, name = "Taxonomic group") +
  scale_x_continuous(trans = pl, breaks = c(0, 1, 3, 10, 30, 100)) +
  scale_y_continuous(trans = pl, breaks = c(0, 1, 3, 10, 30, 100)) +
  facet_wrap(~panel, ncol = 3, labeller = label_wrap_gen(width = 30)) +
  labs(x = "Spread rate", y = expression("Spread rate, second method (km yr"^-1*")")) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", strip.text = element_text(size = 9), legend.text = element_text(size = 9.5), legend.title = element_text(size = 10))
save_fig("figS4_method_agreement", p, 7.4, 4)

# Figure S5
jk <- tb("jackknife.csv")
jk[, method := factor(method, levels = corr)]
full <- cr[predictor == "inat" & response %in% paste0("red_", corr), .(method = factor(sub("red_", "", response), levels = corr), r)]
p <- ggplot(jk, aes(r, dropped)) +
  geom_vline(data = full, aes(xintercept = r), linetype = 2, colour = "#D55E00") +
  geom_vline(xintercept = 0, linetype = 3, colour = "grey60") +
  geom_pointrange(aes(xmin = lwr, xmax = upr), size = 0.25, colour = "#0072B2") +
  facet_wrap(~method, labeller = labeller(method = mlab)) +
  labs(x = expression("Pearson "*italic(r)), y = "Species excluded") +
  theme_bw(base_size = 9) +
  theme(axis.text.y = element_text(face = "italic", size = 8.5), axis.text.x = element_text(size = 9), axis.title = element_text(size = 10),
        strip.text = element_text(size = 9))
save_fig("figS5_jackknife", p, 6.8, 3.8)
