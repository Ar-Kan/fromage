#!/usr/bin/env Rscript

package_version <- function() {
  desc <- read.dcf("DESCRIPTION")
  desc[1, "Version"]
}

run_rcmd <- function(args) {
  cmd <- file.path(R.home("bin"), "Rcmd.exe")
  status <- system2(cmd, args)

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
run_rcmd(c("build", ".", "--no-build-vignettes", "--no-manual"))

message("Checking source package...")
run_rcmd(c("check", tarball, "--no-manual", "--ignore-vignettes"))

message("Full package check completed.")
