#!/usr/bin/env Rscript
# regenerates documentation and runs the tests

message("Documenting package...")
roxygen2::roxygenise(".")

message("Running tests...")
devtools::test()

message("Development checks completed.")
