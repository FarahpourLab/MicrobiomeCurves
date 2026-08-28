skip_if_not_installed("SummarizedExperiment")
skip_if_not_installed("TreeSummarizedExperiment")

test_that("mc_as_demo_se builds a rectangular object", {
    se <- mc_as_demo_se()

    expect_s4_class(se, "TreeSummarizedExperiment")
    expect_equal(nrow(se), 12)
    expect_equal(ncol(se), 70)
    expect_equal(SummarizedExperiment::assayNames(se), "clr")
    expect_true(
        all(c("subject", "timepoint") %in%
            colnames(SummarizedExperiment::colData(se)))
    )
})

test_that("imputing a SummarizedExperiment adds an assay, not replaces one", {
    se <- mc_as_demo_se()
    out <- suppressWarnings(
        mc_run(se, subject_col = "subject", time_col = "timepoint",
            K = 1, verbose = FALSE)
    )

    expect_equal(SummarizedExperiment::assayNames(out), c("clr", "imputed"))

    # the input assay must survive byte-for-byte
    expect_identical(
        SummarizedExperiment::assay(se, "clr"),
        SummarizedExperiment::assay(out, "clr")[, colnames(se), drop = FALSE]
    )

    # and the imputed assay must be complete
    expect_false(anyNA(SummarizedExperiment::assay(out, "imputed")))
})

test_that("the container route agrees with the table route", {
    # Both entry points must give the same values.
    se <- mc_as_demo_se()
    out <- suppressWarnings(
        mc_run(se, subject_col = "subject", time_col = "timepoint",
            K = 1, seed = 123, verbose = FALSE)
    )
    imp <- SummarizedExperiment::assay(out, "imputed")

    tab <- mc_demo_table()
    run <- suppressWarnings(
        run_wide_as_meta(tab, K = 1, seed = 123, verbose = FALSE)
    )

    # The two routes name their columns differently, so align them on the
    # subject and time each column stands for before comparing.
    cd <- SummarizedExperiment::colData(out)
    se_key <- paste(cd$subject, cd$timepoint, sep = "
")

    run_key <- paste(
        run$metadata$subject, run$metadata$time,
        sep = "
"
    )
    cols <- setdiff(names(run$completed), "taxon")
    run_key <- run_key[match(cols, run$metadata$sample)]

    shared <- intersect(se_key, run_key)
    expect_gt(length(shared), 0)

    a <- imp[, match(shared, se_key), drop = FALSE]
    b <- as.matrix(run$completed[, cols[match(shared, run_key)], drop = FALSE])

    expect_equal(unname(a), unname(b), tolerance = 1e-10)
})
