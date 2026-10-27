# Simulated invasions with known radial spread over the recording effort of four target groups.
source("R/functions.R")
a <- as.integer(commandArgs(trailingOnly = TRUE))
part <- if (length(a) > 0) a[1] else 1L
parts <- if (length(a) > 1) a[2] else 1L
dir.create("output/sim", showWarnings = FALSE, recursive = TRUE)

design <- expand.grid(tg = c("Bivalvia", "Perciformes", "Amphipoda", "Cambaridae"), v = c(3, 8),
                      t0 = c(1950L, 1980L), scenario = c("effort", "boost"), rep = 1:4, stringsAsFactors = FALSE)
design$id <- seq_len(nrow(design))
q0 <- 0.02
cache <- list()

for (i in seq(part, nrow(design), by = parts)) {
  z <- design[i, ]
  out <- sprintf("output/sim/sim_%03d.rds", z$id)
  if (file.exists(out)) next
  if (is.null(cache[[z$tg]])) {
    g <- grid(readRDS(taxon_file(z$tg))[year %in% years, .(X, Y, year)], 25)
    dom <- unique(g[, .(cell, gx, gy)])
    dom[, `:=`(cx = (gx + 0.5) * 25, cy = (gy + 0.5) * 25)]
    dom <- merge(dom, g[, .N, by = cell], by = "cell")
    cache[[z$tg]] <- list(g = g, dom = dom)
  }
  g <- copy(cache[[z$tg]]$g)
  dom <- copy(cache[[z$tg]]$dom)

  set.seed(880000L + z$id)
  cand <- dom[N >= median(N)]
  s0 <- cand[sample.int(nrow(cand), 1)]
  dom[, d := sqrt((cx - s0$cx)^2 + (cy - s0$cy)^2)]
  g[, d := dom$d[match(cell, dom$cell)]]
  rt <- z$v * (g$year - z$t0)
  q <- if (z$scenario == "effort") q0 else q0 * (1 + 4 * plogis((g$year - 2005) / 3))
  f <- g[rt >= 0 & d <= rt & runif(nrow(g)) < q, .(X, Y, year)]
  if (nrow(f) < 20) {
    saveRDS(list(design = z, skipped = TRUE, n_focal = nrow(f)), out)
    next
  }
  x <- run(f, g[, .(X, Y, year)], 25, n_rep = 50L)

  truth <- data.table(year = years, aoo = vapply(years, function(y) {
    r <- z$v * (y - z$t0)
    if (r < 0) 0L else sum(dom$d <= r)
  }, integer(1)))
  res <- rates(x$curves, 25)
  res[, true_rate := mapply(function(lo, hi) {
    tr <- radial(truth, 25, c(lo, hi))
    if (is.null(tr)) NA_real_ else slope(tr)
  }, win_lo, win_hi)]
  res[, `:=`(id = z$id, tg = z$tg, v = z$v, t0 = z$t0, scenario = z$scenario, rep = z$rep, n_focal = nrow(f))]
  saveRDS(list(design = z, skipped = FALSE, rates = res, curves = rbindlist(x$curves, idcol = "method"), truth = truth), out)
  message(format(Sys.time(), "%H:%M"), " sim ", z$id)
}
