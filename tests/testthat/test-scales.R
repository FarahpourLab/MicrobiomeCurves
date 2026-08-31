# The completed table is written on three scales. CLR is what the model
# produced. The other two are derived, and how faithful they can be depends
# on what the caller supplied.

make_counts <- function() {
    set.seed(9)
    subs <- paste0("SUB", sprintf("%02d", 1:6))
    days <- c(0, 7, 14)
    meta <- expand.grid(
        SubjectID = subs, Day = days,
        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
    meta <- meta[order(meta$SubjectID, meta$Day), ]
    meta$SampleID <- sprintf("S%03d", seq_len(nrow(meta)))
    meta <- meta[, c("SampleID", "SubjectID", "Day")]
    meta <- meta[!(meta$SubjectID == "SUB03" & meta$Day == 14), ]

    taxa <- paste0("Genus_", LETTERS[1:6])
    counts <- matrix(
        stats::rpois(length(taxa) * nrow(meta), c(400, 200, 90, 40, 15, 5)),
        nrow = length(taxa), dimnames = list(taxa, meta$SampleID)
    )
    counts[5:6, 1:4] <- 0
    list(counts = counts, metadata = meta)
}

run_it <- function(x, type, dir) {
    d <- make_counts()
    suppressWarnings(mc_run(
        x, d$metadata,
        sample_col = "SampleID", subject_col = "SubjectID", time_col = "Day",
        abundance_type = type, out_dir = dir, plots = FALSE,
        C = 1, verbose = FALSE
    ))
}

test_that("inverting a CLR gives proportions that sum to one", {
    x <- c(2, 1, -1, -2)
    ra <- mc_clr_to_ra(x)

    expect_equal(sum(ra), 1)
    expect_true(all(ra > 0))
    # Order is preserved: the largest CLR is the most abundant taxon.
    expect_equal(order(ra), order(x))
})

test_that("a value below the threshold comes back as a zero", {
    x <- c(3, 2, 1, -6)
    plain <- mc_clr_to_ra(x)
    zeroed <- mc_clr_to_ra(x, pseudocount = plain[4] * 2)

    expect_equal(zeroed[4], 0)
    expect_equal(sum(zeroed), 1)
    # The mass taken off is spread over the taxa that remain, so each of
    # them gains.
    expect_true(all(zeroed[1:3] > plain[1:3]))
})

test_that("an all-missing sample stays missing on every scale", {
    x <- rep(NA_real_, 5)

    expect_true(all(is.na(mc_clr_to_ra(x))))
    expect_true(is.na(mc_implied_depth(x)))
    expect_true(all(is.na(mc_ra_to_counts(c(0.5, 0.5), NA_real_))))
})

test_that("implied depth is the reciprocal of the rarest taxon, capped", {
    expect_equal(mc_implied_depth(c(0.5, 0.4, 0.1)), 10L)
    expect_equal(mc_implied_depth(c(0.999999, 0.000001)), 100000L)
    expect_equal(mc_implied_depth(c(0.5, 0.5, 0)), 2L)
})

test_that("counts are the same every time the same data is run", {
    ra <- c(0.5, 0.3, 0.2)

    expect_equal(mc_ra_to_counts(ra, 1000), c(500, 300, 200))
    expect_identical(mc_ra_to_counts(ra, 731), mc_ra_to_counts(ra, 731))
})

test_that("three tables and a log are written", {
    dir <- withr::local_tempdir()
    d <- make_counts()
    run <- run_it(d$counts, "raw", dir)

    expect_setequal(
        basename(run$files),
        c(
            "imputed_clr.tsv", "imputed_relative_abundance.tsv",
            "imputed_counts.tsv", "imputation_log.txt"
        )
    )
    expect_true(all(file.exists(run$files)))
})

test_that("counts in, a collected sample is written back as it was given", {
    dir <- withr::local_tempdir()
    d <- make_counts()
    run <- run_it(d$counts, "raw", dir)

    cn <- utils::read.delim(
        file.path(dir, "imputed_counts.tsv"),
        check.names = FALSE
    )
    ra <- utils::read.delim(
        file.path(dir, "imputed_relative_abundance.tsv"),
        check.names = FALSE
    )
    obs <- colnames(d$counts)
    truth <- d$counts[match(cn$taxon, rownames(d$counts)), ]

    expect_equal(as.matrix(cn[, obs]), truth, ignore_attr = TRUE)
    expect_equal(
        as.matrix(ra[, obs]),
        sweep(truth, 2, colSums(truth), "/"),
        ignore_attr = TRUE
    )
    # Including the structural zeros, which the back-transform alone would
    # have turned into small positive numbers.
    expect_equal(sum(as.matrix(ra[, obs]) == 0), sum(d$counts == 0))
})

test_that("counts in, a created sample sits on the study's scale", {
    dir <- withr::local_tempdir()
    d <- make_counts()
    run <- run_it(d$counts, "raw", dir)

    cn <- utils::read.delim(
        file.path(dir, "imputed_counts.tsv"),
        check.names = FALSE
    )
    made <- setdiff(names(cn), c("taxon", colnames(d$counts)))

    expect_length(made, 1)
    expect_equal(
        sum(cn[[made]]),
        stats::median(colSums(d$counts)),
        tolerance = 1e-6
    )
})

test_that("CLR in, zeros cannot be recovered and the log says so", {
    dir <- withr::local_tempdir()
    d <- make_counts()
    run <- run_it(mc_clr(d$counts, "auto"), "clr", dir)

    ra <- utils::read.delim(
        file.path(dir, "imputed_relative_abundance.tsv"),
        check.names = FALSE
    )
    expect_true(all(as.matrix(ra[, -1]) > 0))
    expect_equal(unname(colSums(ra[, -1])), rep(1, ncol(ra) - 1))

    log <- readLines(file.path(dir, "imputation_log.txt"))
    expect_true(any(grepl("small positive number", log)))
})

test_that("CLR in, one library size is used for the whole table", {
    dir <- withr::local_tempdir()
    d <- make_counts()
    run <- run_it(mc_clr(d$counts, "auto"), "clr", dir)

    cn <- utils::read.delim(
        file.path(dir, "imputed_counts.tsv"),
        check.names = FALSE
    )
    depths <- colSums(cn[, -1])

    # Equal up to the rounding of each column, so samples stay comparable.
    expect_lt(max(depths) - min(depths), nrow(cn))
})

test_that("the CLR table is what the model produced, unaltered", {
    dir <- withr::local_tempdir()
    d <- make_counts()
    run <- run_it(d$counts, "raw", dir)

    written <- utils::read.delim(
        file.path(dir, "imputed_clr.tsv"),
        check.names = FALSE
    )
    expect_equal(written$taxon, run$completed$taxon)
    expect_equal(
        as.matrix(written[, -1]),
        as.matrix(run$completed[, -1]),
        ignore_attr = TRUE, tolerance = 1e-7
    )
})
