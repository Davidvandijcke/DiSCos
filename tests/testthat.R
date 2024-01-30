# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

library(testthat)
library(DiSCos)

options(datatable.auto.thread = 1)  # Limit data.table to 1 thread
Sys.setenv(OMP_NUM_THREADS = 1) # limit environment to 1 thread

test_check("DiSCos")
