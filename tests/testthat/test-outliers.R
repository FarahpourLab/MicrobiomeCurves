# mc_detect_outliers().
#
# The depth is fda.usc's, so these do not re-test it. What is tested is
# that the report is shaped as promised, that the flag means the alpha
# quantile and is reported as a rank rather than a verdict, that depth is
# relative to the subjects in play, and that the screening reported here
# is the screening the imputation applies.

mc_od <- function(n_sub = 12, n_t = 6) {
    d <- expand.grid(
        time = seq_len(n_t), subject = paste0("S", seq_len(n_sub)),
        stringsAsFactors = FALSE
    )
    d <- d[order(d$subject, d$time), ]
    d$sample <- paste0(d$subject, "_T", d$time)
    rownames(d) <- NULL
    d
}

mc_ot <- function(meta, odd = character(0), seed = 3) {
    withr::with_seed(seed, {
        m <- matrix(
            stats::rnorm(2 * nrow(meta), sd = 0.2),
            nrow = 2,
            dimnames = list(c("A", "B"), meta$sample)
        )
        for (nm in rownames(m)) {
            m[nm, ] <- m[nm, ] + sin(meta$time)
        }
        # One subject sent far from the bundle, so there is something a
        # depth ordering ought to find.
        if (length(odd) > 0) {
            hit <- meta$subject %in% odd
            m["A", hit] <- m["A", hit] + 6
        }
        m
    })
}


test_that("every trajectory is reported with its depth", {
    meta <- mc_od()
    tab <- mc_ot(meta)

    od <- mc_detect_outliers(
        tab, meta, "sample", "subject", "time",
        verbose = FALSE
    )

    expect_s3_class(od, "mc_outliers")
    expect_identical(
        names(od$results), c("taxon", "subject", "depth", "outlier")
    )
    expect_identical(nrow(od$results), 2L * 12L)
    expect_type(od$results$depth, "double")
    expect_type(od$results$outlier, "logical")
    expect_identical(
        names(od$n_outliers), c("taxon", "n_subjects", "n_outliers")
    )
})

test_that("the subject sent far from the bundle is the least central", {
    meta <- mc_od()
    tab <- mc_ot(meta, odd = "S7")

    od <- mc_detect_outliers(
        tab, meta, "sample", "subject", "time",
        verbose = FALSE
    )
    a <- od$results[od$results$taxon == "A", ]
    expect_identical(a$subject[which.min(a$depth)], "S7")
    expect_true(a$outlier[a$subject == "S7"])
})

test_that("alpha sets how many are flagged, because it is a quantile", {
    meta <- mc_od(n_sub = 20)
    tab <- mc_ot(meta)

    few <- mc_detect_outliers(tab, meta, "sample", "subject", "time",
        alpha = 0.05, verbose = FALSE
    )
    many <- mc_detect_outliers(tab, meta, "sample", "subject", "time",
        alpha = 0.3, verbose = FALSE
    )

    # The flag is a rank, so a larger alpha flags more of the same curves
    # rather than finding anything new.
    expect_gt(
        sum(many$results$outlier), sum(few$results$outlier)
    )
    flagged_few <- few$results$subject[few$results$outlier]
    flagged_many <- many$results$subject[many$results$outlier]
    expect_true(all(flagged_few %in% flagged_many))
})

test_that("the report is the screening the imputation applies", {
    # Both go through detect_outliers_depth() at the same alpha, so a
    # subject flagged here is the subject the fit would set aside.
    meta <- mc_od()
    tab <- mc_ot(meta, odd = "S4")

    od <- mc_detect_outliers(tab, meta, "sample", "subject", "time",
        taxa = "A", verbose = FALSE
    )

    i <- which(meta$subject == meta$subject)
    Ly <- split(as.numeric(tab["A", ]), meta$subject)
    Lt <- split(meta$time, meta$subject)
    direct <- detect_outliers_depth(Ly, Lt, alpha = 0.05)

    expect_identical(
        od$results$outlier[match(names(Ly), od$results$subject)],
        unname(direct)
    )
    expect_length(i, nrow(meta))
})

test_that("depth is relative to the subjects in play", {
    meta <- mc_od()
    tab <- mc_ot(meta, odd = "S7")

    all_subj <- mc_detect_outliers(tab, meta, "sample", "subject", "time",
        taxa = "A", verbose = FALSE
    )
    without <- mc_detect_outliers(tab, meta, "sample", "subject", "time",
        taxa = "A", subjects = setdiff(paste0("S", 1:12), "S7"),
        verbose = FALSE
    )

    expect_false("S7" %in% without$results$subject)
    # Removing the curve that sat furthest out moves the others' depths:
    # depth is a position among the curves present, not a property of one.
    shared <- intersect(all_subj$results$subject, without$results$subject)
    expect_false(identical(
        all_subj$results$depth[match(shared, all_subj$results$subject)],
        without$results$depth[match(shared, without$results$subject)]
    ))
})

test_that("the report is written when asked for", {
    meta <- mc_od()
    tab <- mc_ot(meta)
    out <- withr::local_tempdir()

    od <- mc_detect_outliers(tab, meta, "sample", "subject", "time",
        out_dir = out, verbose = FALSE
    )
    expect_identical(basename(od$files), "outliers.tsv")
    expect_identical(nrow(utils::read.delim(od$files)), nrow(od$results))
})

test_that("an input the function cannot honour is refused", {
    meta <- mc_od()
    tab <- mc_ot(meta)
    one <- function(...) {
        mc_detect_outliers(tab, meta, "sample", "subject", "time", ...,
            verbose = FALSE
        )
    }

    expect_error(one(alpha = 0), "greater than 0")
    expect_error(one(alpha = 1), "less than 1")
    expect_error(one(alpha = c(0.05, 0.1)), "single number")
    expect_error(one(subjects = "nobody"), "No sample for subject")
})

test_that("printing says the flag is a rank and not a test", {
    meta <- mc_od()
    tab <- mc_ot(meta)
    od <- mc_detect_outliers(tab, meta, "sample", "subject", "time",
        verbose = FALSE
    )

    expect_output(print(od), "mc_outliers")
    expect_output(print(od), "Fraiman-Muniz")
    expect_output(print(od), "not a test")
    expect_identical(withr::with_output_sink(tempfile(), print(od)), od)
})
