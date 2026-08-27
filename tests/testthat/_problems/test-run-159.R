# Extracted from test-run.R:159

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "TaxaTimeImpute", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
run <- quiet_run(taxa_demo, taxon_col = "OTU_ID", K = 1, seed = 42)
dat_full <- taxa_demo
dat_full[["S04.2"]] <- NA_real_
prep <- tti_prepare(
        dat = dat_full,
        taxon_col = "OTU_ID",
        mask_list = run$missing[, c("rep", "time")]
    )
