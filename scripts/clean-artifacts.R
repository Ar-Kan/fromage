#!/usr/bin/env Rscript
# removes only the generated check, archive, and example artifacts listed in that script

paths <- c(
  "fromage.Rcheck",
  Sys.glob("fromage_*.tar.gz"),
  file.path("results", "single-series-example.png"),
  file.path("results", "single-series-example.rds")
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
