suppressPackageStartupMessages({
  library(data.table)
  library(sf)
})

years <- 1950:2024
reps <- 200L
sp <- fread("data/species.csv")
pc <- sp[native_box == 1, species]
meth <- c("naive", "rarefied", "accessible", "listlength")

taxon_file <- function(x) file.path("data", "occ", paste0(gsub(" ", "_", x), ".rds"))

native_box <- function() {
  b <- st_as_sfc(st_bbox(c(xmin = 28, ymin = 40, xmax = 50, ymax = 48), crs = st_crs(4326)))
  st_transform(st_segmentize(b, units::set_units(10000, "m")), 3035)
}

# cleaned records (lon/lat) -> EPSG:3035 with iNaturalist, observation and native-range flags
prep <- function(x, inat = "50c9509d-22c7-4a22-a47d-8c48425ef4a7") {
  xy <- st_coordinates(st_transform(st_as_sf(x, coords = c("lon", "lat"), crs = 4326L), 3035L))
  d <- data.table(X = xy[, 1], Y = xy[, 2], year = x$year, inat = x$datasetKey %in% inat,
                  obs = x$basisOfRecord %in% c("HUMAN_OBSERVATION", "MACHINE_OBSERVATION"))
  p <- st_as_sf(d[, .(X, Y)], coords = c("X", "Y"), crs = 3035L)
  d[, box := lengths(st_intersects(p, native_box())) > 0]
  d[]
}

records <- function(s) {
  f <- readRDS(taxon_file(s))
  g <- readRDS(taxon_file(sp[species == s, target]))
  if (s %in% pc) {
    f <- f[box == FALSE]
    g <- g[box == FALSE]
  }
  list(f = f[, .(X, Y, year)], g = g[, .(X, Y, year)])
}

grid <- function(d, km) {
  d <- copy(d)
  d[, gx := as.integer(floor(X / (km * 1000)))]
  d[, gy := as.integer(floor(Y / (km * 1000)))]
  d[, cell := paste(gx, gy, sep = "_")]
  d[]
}

# cumulative number of cells by year of first record
cumcells <- function(d) cumsum(tabulate(match(d[, min(year), by = cell]$V1, years), length(years)))

quota <- function(f, g) as.integer(quantile(cumsum(table(c(f$year, g$year))), 0.5, names = FALSE))

rarefy <- function(f, g, q, n_rep = reps) {
  p <- rbindlist(list(f[, .(cell, year, foc = TRUE)], g[, .(cell, year, foc = FALSE)]))
  setkey(p, year)
  y <- years[years >= min(f$year)]
  n <- findInterval(y, p$year)
  m <- matrix(NA_real_, length(y), n_rep)
  for (i in which(n >= q)) {
    k <- seq_len(n[i])
    cl <- p$cell[k]
    fo <- p$foc[k]
    m[i, ] <- vapply(seq_len(n_rep), function(b) {
      s <- sample.int(n[i], q, replace = FALSE)
      uniqueN(cl[s][fo[s]])
    }, numeric(1))
  }
  list(year = y, m = m)
}

rcurve <- function(r, mono = TRUE) {
  ok <- !is.na(r$m[, 1])
  a <- apply(r$m[ok, , drop = FALSE], 1, quantile, probs = 0.5, names = FALSE)
  data.table(year = r$year[ok], aoo = if (mono) cummax(a) else a)
}

accessible <- function(f, g, km, radius) {
  m <- km * 1000
  o <- unique(f[, .(gx, gy)])
  ox <- (o$gx + 0.5) * m
  oy <- (o$gy + 0.5) * m
  tc <- unique(g[, .(cell, gx, gy)])
  tx <- (tc$gx + 0.5) * m
  ty <- (tc$gy + 0.5) * m
  dmin <- vapply(seq_along(tx), function(i) min(sqrt((ox - tx[i])^2 + (oy - ty[i])^2)), numeric(1))
  g[cell %in% tc$cell[dmin <= radius * 1000]]
}

listlength <- function(f, g, cumulative = FALSE) {
  d <- merge(g[, .(n = .N), by = .(cell, year)], unique(f[, .(cell, year)])[, p := 1L],
             by = c("cell", "year"), all.x = TRUE)
  d[is.na(p), p := 0L]
  d <- d[year %in% years]
  empty <- data.table(year = integer(), aoo = numeric())
  if (uniqueN(d$year) < 5 || sum(d$p) < 5) return(empty)
  m <- tryCatch(glm(p ~ log1p(n) + factor(year), family = binomial(), data = d), error = function(e) NULL)
  if (is.null(m)) return(empty)
  x <- data.table(year = sort(unique(d$year)), n = median(d$n))
  x[, p := predict(m, newdata = x, type = "response")]
  k <- if (cumulative) cumcells(g)[match(x$year, years)] else uniqueN(d$cell)
  x[, aoo := cummax(p * k)]
  x <- merge(data.table(year = years), x[, .(year, aoo)], by = "year", all.x = TRUE)
  x[, aoo := nafill(aoo, "locf")]
  x[!is.na(aoo)]
}

# two-breakpoint continuous piecewise-linear fit to log(1 + cells), exhaustive search
breaks <- function(d, k = 4L) {
  d <- d[aoo > 0]
  n <- nrow(d)
  if (n < 3 * k) return(c(NA_real_, NA_real_))
  x <- d$year
  y <- log1p(d$aoo)
  cand <- x[k:(n - k + 1)]
  best <- Inf
  bp <- c(NA_real_, NA_real_)
  for (b1 in cand) for (b2 in cand[cand - b1 >= k]) {
    e <- lm.fit(cbind(1, x, pmax(0, x - b1), pmax(0, x - b2)), y)$residuals
    s <- sum(e * e)
    if (s < best - 1e-12) {
      best <- s
      bp <- c(b1, b2)
    }
  }
  bp
}

radial <- function(d, km, w) {
  d <- d[year >= w[1] & year <= w[2] & aoo > 0]
  if (nrow(d) < 3) return(NULL)
  d[, r := sqrt(aoo * km^2 / pi)]
  d[]
}

slope <- function(x) lm.fit(cbind(1, x$year), x$r)$coefficients[[2]]

fit <- function(d, km, ci = TRUE) {
  bp <- breaks(d)
  w <- bp
  if (anyNA(bp)) {
    r <- range(d[aoo > 0, year])
    w <- r + c(0.2, -0.2) * diff(r)
  }
  x <- radial(d, km, w)
  b <- if (is.null(x)) NA_real_ else slope(x)
  if (!ci) return(b)
  lo <- hi <- NA_real_
  if (!is.null(x) && nrow(x) >= 4) {
    m <- lm(r ~ year, data = x)
    v <- tryCatch(suppressWarnings(sandwich::NeweyWest(m, lag = floor(4 * (nrow(x) / 100)^(2 / 9)), prewhite = FALSE, adjust = TRUE)),
                  error = function(e) NULL)
    if (!is.null(v)) {
      e <- qt(0.975, nrow(x) - 2) * sqrt(v[2, 2])
      lo <- b - e
      hi <- b + e
    }
  }
  data.table(rate = b, lwr = lo, upr = hi, lag_end = bp[1], sat_start = bp[2], win_lo = w[1], win_hi = w[2],
             n = if (is.null(x)) 0L else nrow(x), first_year = min(d$year))
}

run <- function(f, g, km, n_rep = reps, radius = 150) {
  f <- grid(f, km)
  g <- grid(g, km)
  q <- quota(f, g)
  rc <- rarefy(f, g, q, n_rep)
  ga <- accessible(f, g, km, radius)
  qa <- if (nrow(ga) > 0) min(q, quota(f, ga)) else q
  ra <- rarefy(f, ga, qa, n_rep)
  cv <- list(naive = data.table(year = years, aoo = cumcells(f)),
             rarefied = rcurve(rc),
             accessible = rcurve(ra),
             listlength = listlength(f, g),
             listlength_cum = listlength(f, g, cumulative = TRUE))
  list(curves = cv, rar = list(rarefied = rc, accessible = ra), quota = c(q, qa))
}

rates <- function(cv, km) rbindlist(lapply(names(cv), function(m) cbind(method = m, fit(cv[[m]], km))))

# rate of each rarefaction replicate (running maximum, breakpoints and regression as for the median curve)
rep_rates <- function(r, km) {
  ok <- !is.na(r$m[, 1])
  apply(r$m[ok, , drop = FALSE], 2, function(v) fit(data.table(year = r$year[ok], aoo = cummax(v)), km, ci = FALSE))
}
