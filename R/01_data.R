# Download the GBIF occurrence downloads listed in data/gbif_downloads.csv, clean and project them.
source("R/functions.R")
library(rgbif)

dl <- fread("data/gbif_downloads.csv", colClasses = "character")
dir.create("data/raw", showWarnings = FALSE, recursive = TRUE)
dir.create("data/occ", showWarnings = FALSE, recursive = TRUE)

basis <- c("HUMAN_OBSERVATION", "MACHINE_OBSERVATION", "PRESERVED_SPECIMEN", "OBSERVATION", "OCCURRENCE")
cols <- c("gbifID", "species", "decimalLatitude", "decimalLongitude", "coordinateUncertaintyInMeters",
          "year", "basisOfRecord", "datasetKey", "countryCode")
europe <- rnaturalearth::ne_countries(scale = "medium", continent = "Europe", returnclass = "sf") |>
  st_make_valid() |>
  st_union()

read_download <- function(key) {
  z <- file.path("data/raw", paste0(key, ".zip"))
  if (!file.exists(z)) occ_download_get(key, path = "data/raw", overwrite = TRUE)
  ex <- file.path(tempdir(), key)
  unzip(z, exdir = ex)
  x <- fread(list.files(ex, "\\.csv$", full.names = TRUE)[1], select = cols, quote = "", showProgress = FALSE)
  unlink(ex, recursive = TRUE)
  setnames(x, c("decimalLatitude", "decimalLongitude"), c("lat", "lon"))
  x
}

clean <- function(x) {
  x <- x[!is.na(lat) & !is.na(lon) & !is.na(year) & basisOfRecord %in% basis]
  x <- x[is.na(coordinateUncertaintyInMeters) | coordinateUncertaintyInMeters <= 10000]
  x <- unique(x, by = c("lat", "lon", "year", "species"))
  p <- st_as_sf(x, coords = c("lon", "lat"), crs = 4326L, remove = FALSE)
  x <- x[lengths(st_intersects(p, europe)) > 0]
  ok <- CoordinateCleaner::clean_coordinates(as.data.frame(x), lon = "lon", lat = "lat", species = "species",
                                             tests = c("capitals", "centroids", "equal", "gbif", "institutions", "zeros"),
                                             value = "flagged", verbose = FALSE)
  x[ok]
}

for (i in seq_len(nrow(dl))) {
  out <- taxon_file(dl$taxon[i])
  if (file.exists(out)) next
  message(dl$taxon[i])
  saveRDS(prep(clean(read_download(dl$key[i]))), out)
}
