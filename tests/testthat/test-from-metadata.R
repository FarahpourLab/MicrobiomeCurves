study <- function(drop = NULL, blank = NULL) {
    subjects <- paste0("M", sprintf("%02d", 1:5))
    days <- c(0, 7, 14)

    meta <- expand.grid(
        subject = subjects, day = days,
        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
    meta$sample <- paste0(meta$subject, "_d", meta$day)
    meta <- meta[order(meta$subject, meta$day), c("sample", "subject", "day")]
    if (!is.null(drop)) meta <- meta[meta$sample != drop, ]
    rownames(meta) <- NULL

    taxa <- paste0("Taxon", sprintf("%02d", 1:4))
    set.seed(42)
    m <- matrix(
        rnorm(length(taxa) * nrow(meta)),
        nrow = length(taxa), dimnames = list(taxa, meta$sample)
    )
    if (!is.null(blank)) m[, blank] <- NA_real_

    list(
        abundance = data.frame(taxon = taxa, m, check.names = FALSE),
        metadata = meta, taxa = taxa, subjects = subjects, days = days
    )
}

build <- function(s, ...) {
    tti_from_metadata(
        s$abundance, s$metadata,
        sample_col = "sample", subject_col = "subject", time_col = "day",
        taxon_col = "taxon", verbose = FALSE, ...
    )
}

test_that("a complete design reports nothing missing", {
    d <- build(study())
    expect_s3_class(d, "tti_design")
    expect_equal(nrow(d$missing), 0)
    expect_equal(length(d$subjects), 5)
    expect_equal(d$times, c(0, 7, 14))
    expect_equal(nrow(d$map), 15)
})

test_that("a sample absent from the metadata is found from the grid", {
    d <- build(study(drop = "M03_d7"))
    expect_equal(nrow(d$missing), 1)
    expect_equal(d$missing$subject, "M03")
    expect_equal(d$missing$time, 7)
    expect_equal(d$missing$reason, "absent_sample")
    expect_true(is.na(d$missing$sample))
})

test_that("a sample whose column holds no data is reported separately", {
    d <- build(study(blank = "M02_d14"))
    expect_equal(nrow(d$missing), 1)
    expect_equal(d$missing$reason, "no_data")
    expect_equal(d$missing$sample, "M02_d14")
})

test_that("both kinds of gap are found together", {
    d <- build(study(drop = "M03_d7", blank = "M02_d14"))
    expect_equal(nrow(d$missing), 2)
    expect_setequal(d$missing$reason, c("absent_sample", "no_data"))
})

test_that("the metadata columns may be called anything", {
    s <- study()
    names(s$metadata) <- c("Run_ID", "cage_animal", "hours_post_gavage")
    d <- tti_from_metadata(
        s$abundance, s$metadata,
        sample_col = "Run_ID", subject_col = "cage_animal",
        time_col = "hours_post_gavage", taxon_col = "taxon", verbose = FALSE
    )
    expect_equal(nrow(d$map), 15)
})

test_that("taxa can come from row names instead of a column", {
    s <- study()
    ab <- s$abundance[, setdiff(names(s$abundance), "taxon")]
    rownames(ab) <- s$taxa
    d <- tti_from_metadata(
        ab, s$metadata,
        sample_col = "sample", subject_col = "subject", time_col = "day",
        verbose = FALSE
    )
    expect_equal(d$table$OTU_ID, s$taxa)
})

test_that("mismatches between the two tables are named", {
    s <- study()
    bad <- s$metadata
    bad$sample <- sub("^M", "X", bad$sample)
    expect_error(
        tti_from_metadata(s$abundance, bad,
            sample_col = "sample", subject_col = "subject",
            time_col = "day", taxon_col = "taxon", verbose = FALSE
        ),
        "No sample name is shared"
    )

    short <- s$metadata[-1, ]
    expect_error(
        tti_from_metadata(s$abundance, short,
            sample_col = "sample", subject_col = "subject",
            time_col = "day", taxon_col = "taxon", verbose = FALSE
        ),
        "not described in the metadata"
    )
})

test_that("a misnamed metadata column says which argument named it", {
    s <- study()
    expect_error(
        tti_from_metadata(s$abundance, s$metadata,
            sample_col = "sample", subject_col = "mouse",
            time_col = "day", taxon_col = "taxon", verbose = FALSE
        ),
        "subject_col"
    )
})

test_that("a subject sampled twice at one time is refused", {
    s <- study()
    s$metadata$day[2] <- s$metadata$day[1]
    expect_error(
        build(s),
        "appear more than once"
    )
})

test_that("a non-numeric time column is refused", {
    s <- study()
    s$metadata$day <- paste0("day", s$metadata$day)
    expect_error(build(s), "not coercible to numeric")
})

test_that("tti_run returns the caller's sample names", {
    s <- study(drop = "M03_d7", blank = "M02_d14")
    run <- suppressWarnings(tti_run(
        s$abundance,
        metadata = s$metadata, sample_col = "sample",
        subject_col = "subject", time_col = "day",
        taxon_col = "taxon", K = 1, verbose = FALSE
    ))

    cols <- setdiff(names(run$completed), "taxon")
    # Every original sample comes back under its own name.
    expect_true(all(s$metadata$sample %in% cols))
    # The full grid is covered: 5 subjects x 3 times.
    expect_equal(length(cols), 15)
    # The created column is named for its subject and time, and flagged.
    expect_true("M03_7" %in% cols)
    expect_true(run$metadata$imputed[run$metadata$sample == "M03_7"])
})

test_that("observed values are not modified", {
    s <- study(drop = "M03_d7")
    run <- suppressWarnings(tti_run(
        s$abundance,
        metadata = s$metadata, sample_col = "sample",
        subject_col = "subject", time_col = "day",
        taxon_col = "taxon", K = 1, verbose = FALSE
    ))
    keep <- s$metadata$sample
    expect_equal(
        unname(as.matrix(run$completed[, keep])),
        unname(as.matrix(s$abundance[, keep]))
    )
})

test_that("completed columns are ordered by subject then time", {
    s <- study(drop = "M03_d7")
    run <- suppressWarnings(tti_run(
        s$abundance,
        metadata = s$metadata, sample_col = "sample",
        subject_col = "subject", time_col = "day",
        taxon_col = "taxon", K = 1, verbose = FALSE
    ))
    cols <- setdiff(names(run$completed), "taxon")

    # The created sample sits between its neighbours in time, rather than
    # being appended after every other subject.
    m03 <- cols[grepl("^M03", cols)]
    expect_equal(m03, c("M03_d0", "M03_7", "M03_d14"))

    # Subjects stay in blocks, in their original order.
    subj <- sub("_.*$", "", cols)
    expect_equal(subj, rep(paste0("M0", 1:5), each = 3))
})

test_that("metadata requires all three column arguments", {
    s <- study()
    expect_error(
        tti_run(s$abundance, metadata = s$metadata, sample_col = "sample"),
        "must all be named"
    )
})

test_that("the design reports itself when verbose", {
    s <- study(drop = "M03_d7")
    expect_message(
        tti_from_metadata(
            s$abundance, s$metadata,
            sample_col = "sample", subject_col = "subject",
            time_col = "day", taxon_col = "taxon", verbose = TRUE
        ),
        "missing"
    )
})
