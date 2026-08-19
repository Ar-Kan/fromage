#!/usr/bin/env Rscript
# documents, tests, builds, and runs `R CMD check --as-cran`

package_version <- function() {
  desc <- read.dcf("DESCRIPTION")
  desc[1, "Version"]
}

run_rcmd <- function(args) {
  executable <- if (.Platform$OS.type == "windows") "R.exe" else "R"
  command <- file.path(R.home("bin"), executable)
  environment <- if (.Platform$OS.type == "windows") {
    c("LC_ALL=C", "LC_CTYPE=C")
  } else {
    character()
  }
  status <- system2(command, c("CMD", args), env = environment)

  if (!identical(status, 0L)) {
    quit(status = status)
  }
}

message("Documenting package...")
roxygen2::roxygenise(".")

message("Running tests...")
devtools::test()

tarball <- sprintf("fromage_%s.tar.gz", package_version())

message("Building source package...")
run_rcmd(c("build", "."))

message("Checking source package...")
run_rcmd(c("check", "--as-cran", tarball))

message("Full package check completed.")
