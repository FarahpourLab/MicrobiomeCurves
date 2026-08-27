# Extracted from test-run.R:23

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "TaxaTimeImpute", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
run <- quiet_run(taxa_demo, taxon_col = "OTU_ID", K = 1)
observed_cols <- setdiff(names(taxa_demo), c("OTU_ID", "S02.1", "S07.4"))
for (cl in observed_cols) {
        expect_identical(
            run$completed[[cl]], taxa_demo[[cl]],
            info = paste("column", cl, "was modified")
        )
    }
expect_identical(run$completed$OTU_ID, taxa_demo$OTU_ID)
