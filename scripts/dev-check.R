#!/usr/bin/env Rscript

message("Documenting package...")
roxygen2::roxygenise(".")

message("Running tests...")
devtools::test()

message("Development checks completed.")
