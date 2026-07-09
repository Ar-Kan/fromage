#!/usr/bin/env Rscript

paths <- c(
  "fromage.Rcheck",
  Sys.glob("fromage_*.tar.gz"),
  file.path("results", "power-plot.png")
)

paths <- paths[file.exists(paths)]

if (length(paths) == 0) {
  message("No generated artifacts found.")
  quit(status = 0)
}

for (path in paths) {
  message("Removing ", path)
  unlink(path, recursive = TRUE, force = TRUE)
}

message("Cleanup completed.")
