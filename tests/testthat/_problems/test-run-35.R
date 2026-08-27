# Extracted from test-run.R:35

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "TaxaTimeImpute", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
run <- quiet_run(taxa_demo, taxon_col = "OTU_ID", K = 1)
expect_equal(nrow(run$completed), nrow(taxa_demo))
expect_equal(
        names(run$completed)[seq_along(names(taxa_demo))],
        names(taxa_demo)
    )
expect_true("S04.2" %in% names(run$completed))
