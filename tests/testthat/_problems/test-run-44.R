# Extracted from test-run.R:44

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "TaxaTimeImpute", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
run <- quiet_run(taxa_demo, taxon_col = "OTU_ID", K = 1)
for (cl in c("S02.1", "S07.4", "S04.2")) {
        expect_false(anyNA(run$completed[[cl]]), info = paste("column", cl))
        expect_type(run$completed[[cl]], "double")
    }
