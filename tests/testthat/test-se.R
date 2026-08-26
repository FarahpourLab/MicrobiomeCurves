skip_if_not_installed("SummarizedExperiment")
skip_if_not_installed("TreeSummarizedExperiment")

test_that("tti_as_demo_se builds a rectangular object", {
    se <- tti_as_demo_se()

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
    se <- tti_as_demo_se()
    out <- suppressWarnings(
        tti_run(se, subject_col = "subject", time_col = "timepoint",
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
    se <- tti_as_demo_se()
    out <- suppressWarnings(
        tti_run(se, subject_col = "subject", time_col = "timepoint",
            K = 1, seed = 123, verbose = FALSE)
    )
    imp <- SummarizedExperiment::assay(out, "imputed")

    tab <- tti_demo_table()
    run <- suppressWarnings(
        tti_run(tab, taxon_col = "OTU_ID", K = 1, seed = 123, verbose = FALSE)
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
})

test_that("run metadata is attached", {
    se <- tti_as_demo_se()
    out <- suppressWarnings(
        tti_run(se, subject_col = "subject", time_col = "timepoint",
            K = 1, verbose = FALSE)
    )
    m <- S4Vectors::metadata(out)$tti_run

    expect_equal(nrow(m$missing), 3)
    expect_equal(m$n_failed, 0)
    expect_equal(m$assay_imputed, "clr")
    expect_type(m$added_samples, "character")
})

test_that("a subject-timepoint with no sample at all is appended", {
    se <- tti_as_demo_se()
    # drop a sample entirely, so it has no column to be written back into
    keep <- colnames(se) != "S05.3"
    se2 <- se[, keep]
    expect_equal(ncol(se2), 69)

    out <- suppressWarnings(
        tti_run(se2, subject_col = "subject", time_col = "timepoint",
            K = 1, verbose = FALSE)
    )

    expect_equal(ncol(out), 70)
    added <- S4Vectors::metadata(out)$tti_run$added_samples
    expect_length(added, 1)

    cd <- SummarizedExperiment::colData(out)
    expect_equal(as.character(cd[added, "subject"]), "S05")
    expect_equal(as.numeric(cd[added, "timepoint"]), 3)
    expect_false(anyNA(SummarizedExperiment::assay(out, "imputed")[, added]))
})

test_that("subjects and times need not be tidy identifiers", {
    # dots in names and non-integer, unevenly spaced times must both work,
    # because the container route recodes them internally.
    set.seed(3)
    n_taxa <- 5
    subs <- c("Pat.01", "Pat.02", "Pat 03", "Pat-04")
    tps <- c(0, 0.5, 2.5, 7)
    grid <- expand.grid(subject = subs, timepoint = tps,
        stringsAsFactors = FALSE)
    mat <- matrix(
        rnorm(n_taxa * nrow(grid)), nrow = n_taxa,
        dimnames = list(paste0("T", seq_len(n_taxa)),
            paste0("smp", seq_len(nrow(grid))))
    )
    mat[, 3] <- NA_real_

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(clr = mat),
        colData = S4Vectors::DataFrame(
            subject = grid$subject,
            timepoint = grid$timepoint,
            row.names = colnames(mat)
        )
    )

    out <- suppressWarnings(
        tti_run(se, subject_col = "subject", time_col = "timepoint",
            K = 1, verbose = FALSE)
    )

    expect_equal(SummarizedExperiment::assayNames(out), c("clr", "imputed"))
    expect_false(anyNA(SummarizedExperiment::assay(out, "imputed")[, 3]))
    # colData is preserved verbatim
    expect_equal(
        as.character(SummarizedExperiment::colData(out)$subject),
        grid$subject
    )
})

test_that("bad container input is rejected clearly", {
    se <- tti_as_demo_se()

    expect_error(
        tti_run(se, subject_col = "nope", time_col = "timepoint"),
        "colData has no column"
    )
    expect_error(
        tti_run(se, subject_col = "subject", time_col = "timepoint",
            assay_name = "missing_assay"),
        "No assay named"
    )
})
