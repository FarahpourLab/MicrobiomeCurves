# Splitting a run on a group column. The promise is that each group is
# imputed as if it had been handed in on its own, so the test that matters
# is the comparison against doing exactly that.

make_arms <- function() {
    set.seed(42)
    subs <- paste0("S", sprintf("%02d", 1:10))
    meta <- expand.grid(
        subject = subs, time = 1:4,
        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
    meta <- meta[order(meta$subject, meta$time), ]
    meta$sample <- paste0(meta$subject, "_T", meta$time)
    meta$arm <- ifelse(
        meta$subject %in% paste0("S", sprintf("%02d", 1:5)),
        "treated", "control"
    )

    taxa <- paste0("Taxon", 1:8)
    m <- matrix(
        stats::rnorm(length(taxa) * nrow(meta)), nrow = length(taxa),
        dimnames = list(taxa, meta$sample)
    )
    list(abundance = scale(m, center = TRUE, scale = FALSE), metadata = meta)
}

drop_one_each <- function(d) {
    gone <- c(
        which(d$metadata$subject == "S02" & d$metadata$time == 3),
        which(d$metadata$subject == "S08" & d$metadata$time == 2)
    )
    list(
        abundance = d$abundance[, -gone, drop = FALSE],
        metadata = d$metadata[-gone, ]
    )
}

run_split <- function(d, ...) {
    suppressWarnings(mc_run(
        d$abundance, d$metadata, sample_col = "sample",
        subject_col = "subject", time_col = "time", abundance_type = "clr",
        group_col = "arm", verbose = FALSE, ...
    ))
}

run_one_arm <- function(d, arm) {
    keep <- d$metadata$sample[d$metadata$arm == arm]
    suppressWarnings(mc_run(
        d$abundance[, colnames(d$abundance) %in% keep, drop = FALSE],
        d$metadata[d$metadata$sample %in% keep, ],
        sample_col = "sample", subject_col = "subject", time_col = "time",
        abundance_type = "clr", out_dir = NULL,
        verbose = FALSE
    ))
}

sorted <- function(p) p[order(p$species, p$subject, p$time), ]

test_that("a split group matches running that group alone", {
    d <- drop_one_each(make_arms())
    split <- run_split(d, out_dir = NULL)

    for (arm in c("treated", "control")) {
        alone <- sorted(run_one_arm(d, arm)$imputed)
        part <- sorted(split$imputed[split$imputed$arm == arm, ])

        expect_equal(nrow(part), nrow(alone))
        expect_equal(part$imputed_value, alone$imputed_value)
        expect_equal(part$lower, alone$lower)
        expect_equal(part$upper, alone$upper)
    }
})

test_that("splitting changes the answer, so it is not a no-op", {
    d <- drop_one_each(make_arms())
    split <- sorted(run_split(d, out_dir = NULL)$imputed)
    pooled <- suppressWarnings(mc_run(
        d$abundance, d$metadata, sample_col = "sample",
        subject_col = "subject", time_col = "time", abundance_type = "clr",
        out_dir = NULL, verbose = FALSE
    ))
    pooled <- sorted(pooled$imputed)

    expect_equal(nrow(split), nrow(pooled))
    expect_false(isTRUE(all.equal(
        split$imputed_value, pooled$imputed_value
    )))
})

test_that("the merged table keeps every sample in its original order", {
    d <- drop_one_each(make_arms())
    split <- run_split(d, out_dir = NULL)
    original <- colnames(d$abundance)

    expect_equal(
        names(split$completed)[seq_along(original) + 1L], original
    )
    # Two samples were absent, so two columns were created for them.
    expect_equal(ncol(split$completed) - 1L, length(original) + 2L)
})

test_that("each group writes its own files and the merge writes its own", {
    d <- drop_one_each(make_arms())
    dir <- withr::local_tempdir()
    run <- run_split(d, out_dir = dir)

    for (arm in c("treated", "control")) {
        expect_true(dir.exists(file.path(dir, arm)))
        expect_true(file.exists(
            file.path(dir, arm, "imputed_clr.tsv")
        ))
    }
    expect_true(file.exists(file.path(dir, "imputed_clr.tsv")))

    # The column is named `group` whatever the caller's column was
    # called, so the file has one shape for every study.
    u <- utils::read.delim(file.path(dir, "imputed_uncertainty.tsv"))
    expect_true("group" %in% names(u))
    expect_false("arm" %in% names(u))
    expect_setequal(unique(u$group), c("treated", "control"))
})

test_that("a subject may not sit in two groups", {
    d <- drop_one_each(make_arms())
    d$metadata$arm[d$metadata$subject == "S03" & d$metadata$time == 1] <-
        "control"

    expect_error(run_split(d, out_dir = NULL), "more than")
})

test_that("the group column has to exist and be complete", {
    d <- drop_one_each(make_arms())
    expect_error(
        suppressWarnings(mc_run(
            d$abundance, d$metadata, sample_col = "sample",
            subject_col = "subject", time_col = "time",
            abundance_type = "clr", group_col = "nope",
            out_dir = NULL, verbose = FALSE
        )),
        "not in the metadata"
    )

    d$metadata$arm[1] <- NA
    expect_error(run_split(d, out_dir = NULL), "is missing for")
})

test_that("a group with nothing to impute is skipped, not fatal", {
    d <- make_arms()
    gone <- which(d$metadata$subject == "S02" & d$metadata$time == 3)
    d <- list(
        abundance = d$abundance[, -gone, drop = FALSE],
        metadata = d$metadata[-gone, ]
    )

    run <- run_split(d, out_dir = NULL)

    expect_named(run$groups, "treated")
    expect_true(any(grepl("not imputed", run$warnings)))
    # The complete arm keeps its samples even though it was not fitted.
    expect_equal(ncol(run$completed) - 1L, ncol(d$abundance) + 1L)
})

test_that("the uncertainty table names subjects the caller would recognise", {
    d <- drop_one_each(make_arms())
    dir <- withr::local_tempdir()
    run_split(d, out_dir = dir)

    u <- utils::read.delim(file.path(dir, "imputed_uncertainty.tsv"))
    expect_setequal(unique(u$subject), c("S02", "S08"))
    expect_setequal(unique(u$timepoint), c(3, 2))
})
