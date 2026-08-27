# Extracted from test-se.R:57

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "TaxaTimeImpute", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
skip_if_not_installed("SummarizedExperiment")
skip_if_not_installed("TreeSummarizedExperiment")

# test -------------------------------------------------------------------------
se <- tti_as_demo_se()
out <- suppressWarnings(
        tti_run(se, subject_col = "subject", time_col = "timepoint",
            K = 1, seed = 123, verbose = FALSE)
    )
imp <- SummarizedExperiment::assay(out, "imputed")
tab <- tti_demo_table()
run <- suppressWarnings(
        run_wide_as_meta(tab, K = 1, seed = 123, verbose = FALSE)
    )
for (i in seq_len(nrow(run$missing))) {
        lbl <- run$missing$col[i]
        expect_equal(
            unname(imp[, lbl]),
            unname(run$completed[[lbl]]),
            tolerance = 1e-8,
            info = lbl
        )
    }
