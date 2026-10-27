# Accessible-area rarefaction with radii of 75 and 300 km (25 km grid).
source("R/functions.R")
a <- as.integer(commandArgs(trailingOnly = TRUE))
part <- if (length(a) > 0) a[1] else 1L
parts <- if (length(a) > 1) a[2] else 1L
dir.create("output/runs", showWarnings = FALSE, recursive = TRUE)

jobs <- CJ(radius = c(75L, 300L), i = seq_len(nrow(sp)), sorted = FALSE)
for (j in seq(part, nrow(jobs), by = parts)) {
  r <- jobs$radius[j]
  i <- jobs$i[j]
  s <- sp$species[i]
  out <- sprintf("output/runs/%s_r%d.rds", gsub(" ", "_", s), r)
  if (file.exists(out)) next
  d <- records(s)
  f <- grid(d$f, 25)
  g <- grid(d$g, 25)
  ga <- accessible(f, g, 25, r)
  q <- quota(f, g)
  if (nrow(ga) > 0) q <- min(q, quota(f, ga))
  set.seed(5100000L + 97L * i + r)
  z <- rarefy(f, ga, q)
  saveRDS(list(species = s, radius = r, quota = q, rate = fit(rcurve(z), 25, ci = FALSE)), out)
  message(format(Sys.time(), "%H:%M"), " ", s, " ", r, " km")
}
