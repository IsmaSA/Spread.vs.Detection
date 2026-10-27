# Range curves and spread rates, 14 species x 3 grid sizes.
# Optional arguments: part n_parts (e.g. Rscript R/02_rates.R 2 4)
source("R/functions.R")
a <- as.integer(commandArgs(trailingOnly = TRUE))
part <- if (length(a) > 0) a[1] else 1L
parts <- if (length(a) > 1) a[2] else 1L
dir.create("output/runs", showWarnings = FALSE, recursive = TRUE)

jobs <- CJ(km = c(25L, 10L, 50L), i = seq_len(nrow(sp)), sorted = FALSE)
for (j in seq(part, nrow(jobs), by = parts)) {
  km <- jobs$km[j]
  i <- jobs$i[j]
  s <- sp$species[i]
  out <- sprintf("output/runs/%s_%d.rds", gsub(" ", "_", s), km)
  if (file.exists(out)) next
  d <- records(s)
  set.seed(100003L + 7919L * i + km)
  x <- run(d$f, d$g, km)
  saveRDS(list(species = s, km = km, n_focal = nrow(d$f), n_target = nrow(d$g), quota = x$quota,
               rates = rates(x$curves, km),
               free = sapply(x$rar, function(z) fit(rcurve(z, mono = FALSE), km, ci = FALSE)),
               reps = lapply(x$rar, rep_rates, km = km),
               curves = rbindlist(x$curves, idcol = "method")), out)
  message(format(Sys.time(), "%H:%M"), " ", s, " ", km, " km")
}
