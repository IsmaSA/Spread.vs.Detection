for (f in c("01_data.R", "02_rates.R", "03_radius.R", "04_simulation.R", "05_tables.R", "06_figures.R")) {
  message("== ", f)
  source(file.path("R", f), local = new.env())
}
writeLines(capture.output(sessionInfo()), "session_info.txt")
