# mc_cluster().
#
# The clustering itself is the machinery the imputation already uses, so
# these do not re-test k-means or FPCA. What is tested is that each of the
# three questions is answered in its own shape, that the result says what
# it clustered, that a set of curves too small or too alike to split comes
# back as one cluster rather than as an error, and that a run is
# reproducible without standing on the caller's random stream.

mc_cd <- function(n_sub = 8, n_t = 6) {
    d <- expand.grid(
        time = seq_len(n_t), subject = paste0("S", seq_len(n_sub)),
        stringsAsFactors = FALSE
    )
    d <- d[order(d$subject, d$time), ]
    d$sample <- paste0(d$subject, "_T", d$time)
    rownames(d) <- NULL
    d
}

mc_ct <- function(meta, effects, seed = 1) {
    withr::with_seed(seed, {
        m <- matrix(
            stats::rnorm(length(effects) * nrow(meta), sd = 0.3),
            nrow = length(effects),
            dimnames = list(names(effects), meta$sample)
        )
        for (nm in names(effects)) m[nm, ] <- m[nm, ] + effects[[nm]](meta)
        m
    })
}

mc_cfx <- function(meta) {
    # Half the subjects rise and half fall, so a clustering that works has
    # something to find.
    up <- function(d) {
        ifelse(d$subject %in% paste0("S", 1:4), 1, -1) * 0.8 * d$time
    }
    list(
        Rising = up,
        Falling = function(d) -up(d),
        Flat = function(d) rep(0, nrow(d))
    )
}


test_that("subjects are clustered within each taxon", {
    meta <- mc_cd()
    tab <- mc_ct(meta, mc_cfx(meta))

    cl <- mc_cluster(tab, meta, "sample", "subject", "time", verbose = FALSE)

    expect_s3_class(cl, "mc_cluster")
    expect_identical(names(cl$results), c("taxon", "subject", "cluster"))
    # Every taxon is asked about every subject.
    expect_identical(nrow(cl$results), 3L * 8L)
    expect_setequal(cl$results$taxon, rownames(tab))
    expect_setequal(unique(cl$results$subject), unique(meta$subject))
    expect_true(all(cl$results$cluster >= 1L))
    expect_identical(cl$settings$by, "subjects")
})

test_that("taxa are clustered within each subject", {
    meta <- mc_cd()
    tab <- mc_ct(meta, mc_cfx(meta))

    cl <- mc_cluster(tab, meta, "sample", "subject", "time",
        by = "taxa", verbose = FALSE
    )

    expect_identical(names(cl$results), c("subject", "taxon", "cluster"))
    # One clustering per subject, over that subject's taxa.
    expect_identical(nrow(cl$results), 8L * 3L)
    expect_identical(nrow(cl$n_clusters), 8L)
})

test_that("taxa are clustered by their mean trajectory", {
    meta <- mc_cd()
    tab <- mc_ct(meta, mc_cfx(meta))

    cl <- mc_cluster(tab, meta, "sample", "subject", "time",
        by = "taxa_mean", verbose = FALSE
    )

    # One clustering over the cohort, so one row per taxon and no unit
    # column: there is only one unit for it to name.
    expect_identical(names(cl$results), c("taxon", "cluster"))
    expect_identical(nrow(cl$results), 3L)
    expect_setequal(cl$results$taxon, rownames(tab))
})

test_that("a fixed C is not exceeded", {
    meta <- mc_cd()
    tab <- mc_ct(meta, mc_cfx(meta))

    cl <- mc_cluster(tab, meta, "sample", "subject", "time",
        C = 2, verbose = FALSE
    )
    expect_true(all(cl$n_clusters$n_clusters <= 2L))
})

test_that("curves too alike to split come back as one cluster", {
    # Identical curves have one distinct score vector between them, and
    # k-means cannot place two centres on one point. The answer is that
    # they do not separate, not that the run fails.
    meta <- mc_cd(n_sub = 4)
    tab <- mc_ct(meta, list(Same = function(d) 0.5 * d$time), seed = 2)
    tab[1, ] <- rep(seq_len(6) * 0.5, 4)

    cl <- mc_cluster(tab, meta, "sample", "subject", "time",
        C = 3, verbose = FALSE
    )
    expect_identical(unique(cl$results$cluster), 1L)
})

test_that("taxa and subjects can be narrowed", {
    meta <- mc_cd()
    tab <- mc_ct(meta, mc_cfx(meta))
    keep <- paste0("S", 1:4)

    cl <- mc_cluster(tab, meta, "sample", "subject", "time",
        taxa = "Rising", subjects = keep, verbose = FALSE
    )
    expect_identical(unique(cl$results$taxon), "Rising")
    expect_setequal(unique(cl$results$subject), keep)
})

test_that("the assignment is written when asked for", {
    meta <- mc_cd()
    tab <- mc_ct(meta, mc_cfx(meta))
    out <- withr::local_tempdir()

    cl <- mc_cluster(tab, meta, "sample", "subject", "time",
        out_dir = out, verbose = FALSE
    )
    expect_identical(basename(cl$files), "clusters.tsv")
    on_disk <- utils::read.delim(cl$files)
    expect_identical(nrow(on_disk), nrow(cl$results))
})

test_that("a run repeats, and leaves the caller's stream alone", {
    meta <- mc_cd()
    tab <- mc_ct(meta, mc_cfx(meta))
    one <- function() {
        mc_cluster(tab, meta, "sample", "subject", "time",
            verbose = FALSE
        )$results$cluster
    }

    set.seed(99)
    before <- stats::runif(1)
    a <- one()
    set.seed(99)
    invisible(stats::runif(1))
    b <- one()
    expect_identical(a, b)

    set.seed(99)
    expect_identical(before, stats::runif(1))
})

test_that("an input the function cannot honour is refused", {
    meta <- mc_cd()
    tab <- mc_ct(meta, mc_cfx(meta))
    one <- function(...) {
        mc_cluster(tab, meta, "sample", "subject", "time", ...,
            verbose = FALSE
        )
    }

    expect_error(one(by = "bogus"), "should be one of")
    expect_error(one(C = 0), "at least 1")
    expect_error(one(C = c(2, 3)), "single number")
    expect_error(one(subjects = "nobody"), "No sample for subject")
})

test_that("printing says what was clustered", {
    meta <- mc_cd()
    tab <- mc_ct(meta, mc_cfx(meta))
    cl <- mc_cluster(tab, meta, "sample", "subject", "time", verbose = FALSE)

    expect_output(print(cl), "mc_cluster")
    expect_output(print(cl), "subjects, within each taxon")
    expect_output(print(cl), "local to each unit")
    expect_identical(withr::with_output_sink(tempfile(), print(cl)), cl)
})
