# Extracted from test-regression.R:73

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "TaxaTimeImpute", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
a <- quiet_run(taxa_demo, taxon_col = "OTU_ID", K = 1, seed = 7)
