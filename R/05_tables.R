# Tables and summary statistics.
source("R/functions.R")
dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)
wr <- function(x, f) fwrite(x, file.path("output/tables", f))
corr <- meth[-1]
f1 <- function(x) ifelse(abs(x) < 0.05, "0.0", sprintf("%.1f", x))

runs <- lapply(list.files("output/runs", "_(10|25|50)[.]rds$", full.names = TRUE), readRDS)
rt <- rbindlist(lapply(runs, function(x) cbind(species = x$species, km = x$km, x$rates)))
ri <- rbindlist(lapply(runs, function(x) data.table(species = x$species, km = x$km, method = names(x$reps),
  rep_lwr = sapply(x$reps, quantile, 0.025, na.rm = TRUE, names = FALSE),
  rep_upr = sapply(x$reps, quantile, 0.975, na.rm = TRUE, names = FALSE))))
rt <- merge(rt, ri, by = c("species", "km", "method"), all.x = TRUE)
rt[, ci_type := fifelse(is.na(rep_lwr), "Newey-West", "rarefaction replicates")]
rt[, `:=`(ci_lwr = fifelse(is.na(rep_lwr), lwr, rep_lwr), ci_upr = fifelse(is.na(rep_upr), upr, rep_upr))]
rt <- merge(sp[, .(species, group, intro_year)], rt, by = "species")
wr(rt[, .(species, group, km, method, rate, ci_lwr, ci_upr, ci_type, lag_end, sat_start, win_lo, win_hi, n, first_year)],
   "spread_rates.csv")
m25 <- rt[km == 25]
cv25 <- rbindlist(lapply(Filter(function(x) x$km == 25, runs), function(x) cbind(species = x$species, x$curves)))
wr(cv25, "range_curves_25km.csv")

# records, iNaturalist and human/machine-observation shares
occ <- rbindlist(lapply(runs, function(x) data.table(species = x$species, km = x$km, n = x$n_focal)))[km == 25, .(species, n)]
wr(occ[match(sp$species, species)], "occurrences.csv")
cs <- rbindlist(lapply(sp$species, function(s) {
  f <- readRDS(taxon_file(s))
  if (s %in% pc) f <- f[box == FALSE]
  f[, .(n = .N, n_inat = sum(inat), obs = mean(obs)), keyby = year][, species := s]
}))
cs[, inat := n_inat / n]
wr(cs[, .(species, year, n, n_inat, inat, obs)], "inat_share_by_year.csv")
sh <- cs[year >= 2010, .(inat = mean(inat), obs = mean(obs)), by = species]
pool <- cs[, .(inat = sum(n_inat) / sum(n)), by = .(period = fifelse(year < 2010, "1950-2009", as.character(year)))]
wr(pool, "inat_share_pooled.csv")

eff <- rbindlist(lapply(unique(sp$target), function(t)
  data.table(target = t, year = years, cells = cumcells(grid(readRDS(taxon_file(t)), 25)))))
wr(eff, "sampled_cells.csv")

# Table S2
ci_txt <- function(r, l, u, m) ifelse(m == "listlength" & abs(r) < 0.05, "flat (no expansion phase)",
                                      sprintf("%.1f (%.1f to %.1f)", r, l, u))
s2 <- dcast(m25[method %in% meth, .(species, method, txt = ci_txt(rate, ci_lwr, ci_upr, method))], species ~ method, value.var = "txt")
ph <- dcast(m25[method %in% c("naive", "rarefied")], species ~ method, value.var = c("lag_end", "sat_start"))
s2 <- merge(s2[, c("species", meth), with = FALSE], ph, by = "species")
s2[, species := ifelse(sp$intro_year[match(species, sp$species)] < 1950, paste(species, "†"), species)]
wr(s2, "tableS2_rates_phases.csv")

# medians, paired tests, interval overlap
set.seed(20261020)
md <- m25[method %in% c(meth, "listlength_cum"), {
  x <- rate[is.finite(rate)]
  b <- quantile(replicate(10000, median(sample(x, replace = TRUE))), c(0.025, 0.975), names = FALSE)
  .(n = length(x), median = median(x), lwr = b[1], upr = b[2])
}, keyby = method]
wr(md, "median_rates.csv")

w <- dcast(m25, species + group + intro_year ~ method, value.var = "rate")
wx <- rbindlist(lapply(c(corr, "listlength_cum"), function(m) {
  t <- suppressWarnings(wilcox.test(w$naive, w[[m]], paired = TRUE, conf.int = TRUE, exact = FALSE))
  data.table(method = m, n = nrow(w), naive_higher = sum(w$naive > w[[m]]), hl = t$estimate[[1]],
             hl_lwr = t$conf.int[1], hl_upr = t$conf.int[2], p = t$p.value)
}))
wr(wx, "wilcoxon.csv")

ov <- merge(m25[method == "naive", .(species, nl = ci_lwr, nu = ci_upr)],
            m25[method %in% corr, .(species, method, cl = ci_lwr, cu = ci_upr)], by = "species")
ov[, overlap := !(nl > cu | cl > nu)]
wr(ov[, .(n = .N, no_overlap = sum(!overlap)), by = method], "interval_overlap.csv")

pp <- rbindlist(lapply(c(TRUE, FALSE), function(b) {
  x <- w[(intro_year < 1950) == b]
  data.table(before_1950 = b, n = nrow(x), method = meth, median = sapply(meth, function(m) median(x[[m]])),
             naive_higher = sapply(meth, function(m) sum(x$naive > x[[m]])))
}))
wr(pp, "intro_before_after_1950.csv")

# invasion phases: naive curve truncated to the years of the rarefied curve
tr <- rbindlist(lapply(sp$species, function(s) {
  f0 <- cv25[species == s & method == "rarefied", min(year)]
  bp <- breaks(cv25[species == s & method == "naive" & year >= f0])
  data.table(species = s, first_rarefied = f0, lag_naive_trunc = bp[1], sat_naive_trunc = bp[2])
}))
ps <- merge(dcast(m25[method %in% meth], species ~ method, value.var = c("lag_end", "sat_start")), tr, by = "species")
ps[, `:=`(lag_shift = lag_end_rarefied - lag_end_naive, lag_shift_matched = lag_end_rarefied - lag_naive_trunc,
          sat_shift = sat_start_rarefied - sat_start_naive, sat_shift_matched = sat_start_rarefied - sat_naive_trunc,
          lag_after_first = lag_end_rarefied - first_rarefied)]
wr(ps, "phases.csv")

# agreement among corrections
ag <- rbindlist(combn(corr, 2, function(p) {
  t <- suppressWarnings(cor.test(w[[p[1]]], w[[p[2]]], method = "spearman"))
  data.table(a = p[1], b = p[2], rho = t$estimate[[1]], p = t$p.value)
}, simplify = FALSE))
rk <- apply(as.matrix(w[, ..corr]), 2, rank)
rs <- rowSums(rk)
ag[, kendall_w := 12 * sum((rs - mean(rs))^2) / (ncol(rk)^2 * (nrow(rk)^3 - nrow(rk)))]
wr(ag, "method_agreement.csv")

# relative reduction vs iNaturalist and human/machine-observation shares
b <- merge(w, sh, by = "species")
for (m in corr) b[[paste0("red_", m)]] <- (b$naive - b[[m]]) / b$naive
wr(b, "reduction_and_shares.csv")
resp <- c("naive", corr, paste0("red_", corr))
cr <- rbindlist(lapply(c("inat", "obs"), function(p) rbindlist(lapply(resp, function(r) {
  t <- cor.test(b[[p]], b[[r]])
  t2 <- suppressWarnings(cor.test(b[[p]], b[[r]], method = "spearman"))
  data.table(predictor = p, response = r, n = nrow(b), r = t$estimate[[1]], lwr = t$conf.int[1], upr = t$conf.int[2],
             p = t$p.value, rho = t2$estimate[[1]], p_rho = t2$p.value)
}))))
wr(cr, "correlations.csv")
jk <- rbindlist(lapply(corr, function(m) rbindlist(lapply(b$species, function(s) {
  x <- b[species != s]
  t <- cor.test(x$inat, x[[paste0("red_", m)]])
  data.table(method = m, dropped = s, r = t$estimate[[1]], lwr = t$conf.int[1], upr = t$conf.int[2], p = t$p.value)
}))))
wr(jk, "jackknife.csv")

# Table S4: grid size
s4 <- dcast(rt[method %in% meth], species ~ method + km, value.var = "rate")
s4 <- s4[, c("species", paste(rep(meth, each = 3), c(10, 25, 50), sep = "_")), with = FALSE]
wr(s4[, lapply(.SD, function(x) if (is.numeric(x)) f1(x) else x)], "tableS4_grid.csv")

# Table S5: sensitivity analyses
nh <- function(a, b) sprintf("%d of %d", sum(a > b), length(a))
s5 <- list()
for (k in c(10L, 50L)) for (m in meth) {
  a <- rt[km == 25 & method == m][order(species)]
  x <- rt[km == k & method == m][order(species)]
  nk <- rt[km == k & method == "naive"][order(species)]
  s5[[length(s5) + 1]] <- data.table(analysis = paste0("Grid ", k, " km"), method = m, n = nrow(x),
    median_main = median(a$rate), median_sens = median(x$rate), rho = cor(a$rate, x$rate, method = "spearman"),
    naive_higher = if (m == "naive") NA_character_ else nh(nk$rate, x$rate))
}
rad <- rbindlist(lapply(list.files("output/runs", "_r[0-9]+[.]rds$", full.names = TRUE), function(f) {
  x <- readRDS(f)
  data.table(species = x$species, radius = x$radius, rate = x$rate)
}))
wr(rad, "accessible_radius.csv")
for (r in c(75L, 300L)) {
  x <- merge(w[, .(species, naive, accessible)], rad[radius == r], by = "species")
  s5[[length(s5) + 1]] <- data.table(analysis = paste0("Accessible-area radius ", r, " km"), method = "accessible", n = nrow(x),
    median_main = median(x$accessible), median_sens = median(x$rate), rho = cor(x$accessible, x$rate, method = "spearman"),
    naive_higher = nh(x$naive, x$rate))
}
fr <- rbindlist(lapply(runs, function(x) data.table(species = x$species, km = x$km, method = names(x$free), free = x$free)))
fr <- merge(fr, rt[, .(species, km, method, rate)], by = c("species", "km", "method"))
wr(fr, "unconstrained_rarefied.csv")
for (m in c("rarefied", "accessible")) {
  x <- fr[km == 25 & method == m]
  s5[[length(s5) + 1]] <- data.table(analysis = "Rarefied curve not forced to be non-decreasing", method = m, n = nrow(x),
    median_main = median(x$rate), median_sens = median(x$free), rho = cor(x$rate, x$free, method = "spearman"),
    naive_higher = NA_character_)
}
s5[[length(s5) + 1]] <- data.table(analysis = "List-length scaled by cumulative sampled cells", method = "listlength", n = nrow(w),
  median_main = median(w$listlength), median_sens = median(w$listlength_cum), rho = NA_real_, naive_higher = NA_character_)
s5 <- rbindlist(s5)
s5[, `:=`(median_main = sprintf("%.1f", median_main), median_sens = sprintf("%.1f", median_sens),
          rho = ifelse(is.na(rho), "–", sprintf("%.2f", rho)), naive_higher = fcoalesce(naive_higher, "–"))]
for (bb in c(TRUE, FALSE)) {
  x <- w[(intro_year < 1950) == bb]
  s5 <- rbind(s5, data.table(analysis = if (bb) "Introduced before 1950" else "Introduced after 1950",
    method = "naive / rarefied", n = nrow(x), median_main = "–",
    median_sens = sprintf("%.1f / %.1f", median(x$naive), median(x$rarefied)), rho = "–",
    naive_higher = nh(x$naive, x$rarefied)))
}
wr(s5, "tableS5_sensitivity.csv")
wr(data.table(n_rates = nrow(fr), n_negative = sum(fr$free < 0), species_negative = paste(unique(fr[free < 0, species]), collapse = "; ")),
   "unconstrained_negative.csv")

# simulations
sims <- lapply(list.files("output/sim", "[.]rds$", full.names = TRUE), readRDS)
ss <- rbindlist(lapply(Filter(function(x) !x$skipped, sims), `[[`, "rates"))
ss[, ratio := rate / true_rate]
wr(ss, "simulation_runs.csv")
s3 <- ss[is.finite(ratio), .(n = .N, median = median(ratio), q25 = quantile(ratio, 0.25), q75 = quantile(ratio, 0.75),
                             within25 = mean(abs((rate - true_rate) / true_rate) <= 0.25),
                             over25 = mean(ratio > 1.25), under25 = mean(ratio < 0.75)), by = .(scenario, method)]
s3 <- merge(s3, ss[is.finite(rate), .(rho = cor(rate, v, method = "spearman")), by = .(scenario, method)], by = c("scenario", "method"))
s3 <- s3[order(-(scenario == "effort"), match(method, c(meth, "listlength_cum")))]
wr(s3, "tableS3_simulation.csv")
x <- ss[is.finite(ratio) & method %in% meth]
wr(data.table(runs = length(sims), skipped = sum(sapply(sims, `[[`, "skipped")), ratios = nrow(x),
              outside_axis = sum(x$ratio < 0.03 | x$ratio > 25)), "simulation_counts.csv")
