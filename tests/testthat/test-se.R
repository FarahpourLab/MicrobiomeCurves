skip_if_not_installed("SummarizedExperiment")
skip_if_not_installed("TreeSummarizedExperiment")

test_that("mc_as_demo_se builds a rectangular object", {
    se <- demo_se_small()

    expect_s4_class(se, "TreeSummarizedExperiment")
    expect_equal(nrow(se), nrow(taxa_demo))

    # Rectangular is the claim: every subject crossed with every time point,
    # including the samples that do not exist, which come back as all-NA
    # columns rather than being left out.
    cd <- SummarizedExperiment::colData(se)
    expect_equal(
        ncol(se), length(unique(cd$subject)) * length(unique(cd$timepoint))
    )
    expect_equal(SummarizedExperiment::assayNames(se), "clr")
    expect_true(
        all(c("subject", "timepoint") %in%
            colnames(SummarizedExperiment::colData(se)))
    )
})

test_that("imputing a SummarizedExperiment adds an assay, not replaces one", {
    se <- demo_se_small()
    out <- suppressWarnings(
        mc_run(se, subject_col = "subject", time_col = "timepoint",
            C = 1, verbose = FALSE)
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
    se <- demo_se_small()
    out <- suppressWarnings(
        mc_run(se, subject_col = "subject", time_col = "timepoint",
            C = 1, seed = 123, verbose = FALSE)
    )
    imp <- SummarizedExperiment::assay(out, "imputed")

    # The same corner of the cohort the SE fixture was cut to, so the two
    # routes are compared on identical data. mc_demo_table() would load the
    # whole bundled table and the comparison would then be against a
    # different, much larger input.
    cd0 <- SummarizedExperiment::colData(se)
    subs <- unique(cd0$subject)
    cols <- setdiff(names(taxa_demo), "OTU_ID")
    tab <- taxa_demo[, c(
        "OTU_ID", cols[sub("\\..*$", "", cols) %in% subs]
    ), drop = FALSE]

    run <- suppressWarnings(
        run_wide_as_meta(tab, C = 1, seed = 123, verbose = FALSE)
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
