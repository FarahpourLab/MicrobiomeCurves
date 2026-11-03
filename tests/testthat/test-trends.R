# mc_trends() and mc_compare().
#
# The model is mgcv's, so these do not re-test GAM fitting. What is tested is
# that the right model is built, that the two parameterisations behave as
# documented, that an input the function cannot honour is refused rather than
# guessed at, and that the result tables are shaped as promised.

mc_td <- function(n_sub = 12, n_t = 6, arms = 2) {
    d <- expand.grid(
        time = seq_len(n_t), subject = paste0("S", seq_len(n_sub)),
        stringsAsFactors = FALSE
    )
    d <- d[order(d$subject, d$time), ]
    d$sample <- paste0(d$subject, "_T", d$time)
    per <- ceiling(n_sub / arms)
    idx <- ceiling(as.integer(sub("^S", "", d$subject)) / per)
    d$arm <- paste0("arm", LETTERS[pmin(idx, arms)])
    d$sex <- ifelse(
        as.integer(sub("^S", "", d$subject)) %% 2 == 0, "F", "M"
    )
    rownames(d) <- NULL
    d
}

mc_tt <- function(meta, effects, seed = 1) {
    withr::with_seed(seed, {
        m <- matrix(
            stats::rnorm(length(effects) * nrow(meta), sd = 0.3),
            nrow = length(effects),
            dimnames = list(names(effects), meta$sample)
        )
        for (nm in names(effects)) m[nm, ] <- m[nm, ] + effects[[nm]](meta)
        data.frame(taxon = rownames(m), m, check.names = FALSE)
    })
}

test_that("a taxon that rises is found and a flat one is not", {
    meta <- mc_td()
    tab <- mc_tt(meta, list(
        Rising = function(d) 0.8 * d$time,
        Flat = function(d) 0
    ))

    tr <- mc_trends(tab, meta, "sample", "subject", "time", verbose = FALSE)

    expect_s3_class(tr, "mc_trends")
    expect_equal(nrow(tr$results), 2L)
    expect_true(all(tr$results$converged))

    p <- stats::setNames(tr$results$p_value, tr$results$taxon)
    expect_lt(p[["Rising"]], 0.01)
    expect_gt(p[["Flat"]], 0.05)
})

test_that("mc_compare finds a taxon whose arms diverge", {
    meta <- mc_td()
    tab <- mc_tt(meta, list(
        Diverges = function(d) ifelse(d$arm == "armB", 0.9 * d$time, 0),
        Parallel = function(d) ifelse(d$arm == "armB", 2, 0)
    ))

    cm <- mc_compare(tab, meta, "sample", "subject", "time",
        group_col = "arm", verbose = FALSE
    )

    expect_s3_class(cm, "mc_compare")
    p <- stats::setNames(cm$results$p_value, cm$results$taxon)
    # Parallel differs in level, not in shape, so the shape test should not
    # fire on it. Diverges differs in shape and should.
    expect_lt(p[["Diverges"]], 0.01)
    expect_gt(p[["Parallel"]], 0.05)
})

test_that("two groups use the ordered parameterisation and name a level", {
    meta <- mc_td(arms = 2)
    tab <- mc_tt(meta, list(X = function(d) 0.5 * d$time))

    cm <- mc_compare(tab, meta, "sample", "subject", "time",
        group_col = "arm", verbose = FALSE
    )
    expect_identical(cm$settings$group, "ordered")
    expect_identical(cm$settings$reference, "armA")
    expect_identical(cm$results$reference[1], "armA")
    expect_identical(cm$results$n_groups[1], 2L)
})

test_that("three groups default to the global parameterisation", {
    meta <- mc_td(n_sub = 18, arms = 3)
    tab <- mc_tt(meta, list(X = function(d) 0.5 * d$time))

    cm <- mc_compare(tab, meta, "sample", "subject", "time",
        group_col = "arm", verbose = FALSE
    )
    expect_identical(cm$settings$group, "global")
    expect_identical(cm$settings$reference, "global")
    expect_identical(cm$results$n_groups[1], 3L)
})

test_that("reference picks the parameterisation and is validated", {
    meta <- mc_td(n_sub = 18, arms = 3)
    tab <- mc_tt(meta, list(X = function(d) 0.5 * d$time))

    cm <- mc_compare(tab, meta, "sample", "subject", "time",
        group_col = "arm", reference = "armC", verbose = FALSE
    )
    expect_identical(cm$settings$group, "ordered")
    expect_identical(cm$settings$reference, "armC")
    expect_identical(cm$settings$levels[1], "armC")

    expect_error(
        mc_compare(tab, meta, "sample", "subject", "time",
            group_col = "arm", reference = "armZ", verbose = FALSE
        ),
        "available"
    )
})

test_that("with 3 groups the test is the omnibus, whatever the reference", {
    # Three arms, and the one that diverges is NOT the first level after the
    # reference. The ordered parameterisation then holds two difference
    # smooths, and reading only the first would miss this taxon entirely.
    meta <- mc_td(n_sub = 18, arms = 3)
    meta$arm <- c("control", "dose_low", "dose_high")[
        match(meta$arm, c("armA", "armB", "armC"))
    ]
    tab <- mc_tt(meta, list(
        LowDiverges = function(d) ifelse(d$arm == "dose_low", 1.2 * d$time, 0),
        AllSame = function(d) 0.4 * d$time
    ))

    # levels sort to control, dose_high, dose_low, so dose_high is the first
    # difference smooth and dose_low the second.
    ref <- mc_compare(tab, meta, "sample", "subject", "time",
        group_col = "arm", reference = "control", verbose = FALSE
    )
    glo <- mc_compare(tab, meta, "sample", "subject", "time",
        group_col = "arm", reference = "global", verbose = FALSE
    )

    for (cm in list(ref, glo)) {
        p <- stats::setNames(cm$results$p_value, cm$results$taxon)
        expect_lt(p[["LowDiverges"]], 0.01)
        expect_gt(p[["AllSame"]], 0.05)
    }

    # Both parameterisations answer the same omnibus question, so they should
    # agree on which taxa are flagged even though the fits differ.
    expect_identical(
        ref$results$taxon[ref$results$q_value < 0.05],
        glo$results$taxon[glo$results$q_value < 0.05]
    )
})

test_that("the log says which test was used", {
    meta3 <- mc_td(n_sub = 18, arms = 3)
    meta2 <- mc_td(arms = 2)
    tab3 <- mc_tt(meta3, list(X = function(d) 0.5 * d$time))
    tab2 <- mc_tt(meta2, list(X = function(d) 0.5 * d$time))
    dir <- withr::local_tempdir()

    two <- mc_compare(tab2, meta2, "sample", "subject", "time",
        group_col = "arm", verbose = FALSE
    )
    three <- mc_compare(tab3, meta3, "sample", "subject", "time",
        group_col = "arm", reference = "armA", verbose = FALSE
    )

    expect_match(two$settings$test, "difference smooth", fixed = TRUE)
    # Three groups take the omnibus even with a named reference.
    expect_match(three$settings$test, "omnibus", fixed = TRUE)
    expect_identical(three$settings$reference, "armA")
})

test_that("mc_compare requires group_col", {
    meta <- mc_td()
    tab <- mc_tt(meta, list(X = function(d) 0.5 * d$time))

    expect_error(
        mc_compare(tab, meta, "sample", "subject", "time", verbose = FALSE),
        "group_col is required"
    )
    expect_error(
        mc_compare(tab, meta, "sample", "subject", "time",
            group_col = "nope", verbose = FALSE
        ),
        "no column"
    )
})

test_that("a non-numeric time column is refused, not coerced", {
    meta <- mc_td()
    meta$label <- paste0("week_", meta$time)
    tab <- mc_tt(meta, list(X = function(d) 0.5 * d$time))

    expect_error(
        mc_trends(tab, meta, "sample", "subject", "label", verbose = FALSE),
        "not numeric"
    )
})

test_that("k is read from the design and bounded by it", {
    meta <- mc_td(n_t = 5)
    tab <- mc_tt(meta, list(X = function(d) 0.5 * d$time))

    tr <- mc_trends(tab, meta, "sample", "subject", "time", verbose = FALSE)
    # 5 distinct times, so min(5 - 1, 10).
    expect_identical(tr$settings$k, 4L)

    expect_error(
        mc_trends(tab, meta, "sample", "subject", "time",
            k = 9, verbose = FALSE
        ),
        "exceeds"
    )
    expect_error(
        mc_trends(tab, meta, "sample", "subject", "time",
            k = 2, verbose = FALSE
        ),
        "at least 3"
    )
})

test_that("NA values are dropped rather than failing the fit", {
    meta <- mc_td()
    tab <- mc_tt(meta, list(X = function(d) 0.8 * d$time))
    gap <- tab
    gap[1, c("S3_T2", "S7_T4", "S11_T5")] <- NA

    full <- mc_trends(tab, meta, "sample", "subject", "time", verbose = FALSE)
    part <- mc_trends(gap, meta, "sample", "subject", "time", verbose = FALSE)

    expect_true(part$results$converged[1])
    expect_identical(part$results$n_obs[1], full$results$n_obs[1] - 3L)
    # The same taxon, three cells fewer: the answer should move a little and
    # not change in kind.
    expect_lt(part$results$p_value[1], 0.01)
})

test_that("each subject_effect builds the documented terms", {
    expect_identical(
        mc_trend_terms(5, "intercept"),
        c("s(time, k = 5)", "s(subject, bs = \"re\")")
    )
    expect_identical(
        mc_trend_terms(5, "slope"),
        c(
            "s(time, k = 5)", "s(subject, bs = \"re\")",
            "s(subject, time, bs = \"re\")"
        )
    )
    expect_identical(
        mc_trend_terms(5, "smooth"),
        c("s(time, k = 5)", "s(time, subject, bs = \"fs\", k = 5, m = 1)")
    )
    expect_identical(mc_trend_terms(5, "none"), "s(time, k = 5)")

    # The null model for the nested comparison drops only the group smooth.
    expect_false(any(grepl(
        "by = mc_grp",
        mc_trend_terms(5, "intercept", group = "global", group_smooth = FALSE)
    )))
    expect_true(any(grepl(
        "by = mc_grp",
        mc_trend_terms(5, "intercept", group = "global", group_smooth = TRUE)
    )))
})

test_that("every subject_effect runs and covariates are accepted", {
    meta <- mc_td()
    tab <- mc_tt(meta, list(X = function(d) 0.8 * d$time))

    for (se in c("intercept", "slope", "smooth", "none")) {
        tr <- mc_trends(tab, meta, "sample", "subject", "time",
            subject_effect = se, verbose = FALSE
        )
        expect_true(tr$results$converged[1], label = se)
        expect_identical(tr$settings$subject_effect, se)
    }

    tr <- mc_trends(tab, meta, "sample", "subject", "time",
        covariates = "sex", verbose = FALSE
    )
    expect_true(tr$results$converged[1])
    expect_identical(tr$settings$covariates, "sex")
})

test_that("a covariate using a reserved name is refused", {
    meta <- mc_td()
    meta$value <- 1
    tab <- mc_tt(meta, list(X = function(d) 0.5 * d$time))

    expect_error(
        mc_trends(tab, meta, "sample", "subject", "time",
            covariates = "value", verbose = FALSE
        ),
        "used"
    )
})

test_that("smooth subject curves warn when group cannot vary within subject", {
    meta <- mc_td()
    tab <- mc_tt(meta, list(X = function(d) 0.5 * d$time))

    expect_warning(
        mc_compare(tab, meta, "sample", "subject", "time",
            group_col = "arm", subject_effect = "smooth", verbose = FALSE
        ),
        "absorb"
    )
})

test_that("only the shared samples are used", {
    meta <- mc_td()
    tab <- mc_tt(meta, list(X = function(d) 0.8 * d$time))
    tab$extra_sample <- stats::rnorm(nrow(tab))
    meta2 <- meta[meta$sample != "S1_T1", ]

    tr <- mc_trends(tab, meta2, "sample", "subject", "time", verbose = FALSE)
    expect_identical(tr$results$n_obs[1], nrow(meta2))
    expect_identical(nrow(tr$design), nrow(meta2))

    expect_error(
        mc_trends(tab, transform(meta, sample = paste0("z", meta$sample)),
            "sample", "subject", "time", verbose = FALSE
        ),
        "No sample name is in both"
    )
})

test_that("a failed taxon gets a row rather than stopping the run", {
    meta <- mc_td()
    tab <- mc_tt(meta, list(
        Good = function(d) 0.8 * d$time,
        Constant = function(d) 0
    ))
    tab[tab$taxon == "Constant", -1] <- 7

    tr <- mc_trends(tab, meta, "sample", "subject", "time", verbose = FALSE)
    expect_identical(nrow(tr$results), 2L)
    row <- tr$results[tr$results$taxon == "Constant", ]
    expect_false(row$converged)
    expect_true(is.na(row$p_value))
    expect_identical(row$note, "constant")
})

test_that("adjustment is over the taxa that were actually tested", {
    meta <- mc_td()
    tab <- mc_tt(meta, list(
        A = function(d) 0.8 * d$time,
        B = function(d) 0.6 * d$time,
        Bad = function(d) 0
    ))
    tab[tab$taxon == "Bad", -1] <- 1

    tr <- mc_trends(tab, meta, "sample", "subject", "time", verbose = FALSE)
    ok <- tr$results[tr$results$converged, ]
    expect_identical(nrow(ok), 2L)
    # Two tests were carried out, so BH multiplies the smaller p by 2, not 3.
    expect_equal(
        sort(ok$q_value),
        sort(stats::p.adjust(ok$p_value, "BH")),
        tolerance = 1e-12
    )
    expect_true(is.na(tr$results$q_value[tr$results$taxon == "Bad"]))
})

test_that("results carry the promised columns in the promised order", {
    meta <- mc_td()
    tab <- mc_tt(meta, list(X = function(d) 0.8 * d$time))

    tr <- mc_trends(tab, meta, "sample", "subject", "time", verbose = FALSE)
    expect_identical(
        names(tr$results),
        c(
            "taxon", "edf", "F", "p_value", "q_value", "k",
            "n_obs", "n_subjects", "converged", "note"
        )
    )

    cm <- mc_compare(tab, meta, "sample", "subject", "time",
        group_col = "arm", verbose = FALSE
    )
    expect_identical(
        names(cm$results),
        c(
            "taxon", "edf", "F", "p_value", "q_value", "k",
            "n_obs", "n_subjects", "converged", "note",
            "n_groups", "reference"
        )
    )
})

test_that("results are ordered by adjusted p-value", {
    meta <- mc_td()
    tab <- mc_tt(meta, list(
        Weak = function(d) 0.15 * d$time,
        Strong = function(d) 1.5 * d$time,
        Flat = function(d) 0
    ))

    tr <- mc_trends(tab, meta, "sample", "subject", "time", verbose = FALSE)
    expect_identical(tr$results$taxon[1], "Strong")
    q <- tr$results$q_value[!is.na(tr$results$q_value)]
    expect_identical(q, sort(q))
})

test_that("taxa selects a subset and rejects a name that is not there", {
    meta <- mc_td()
    tab <- mc_tt(meta, list(
        A = function(d) 0.8 * d$time,
        B = function(d) 0.5 * d$time
    ))

    tr <- mc_trends(tab, meta, "sample", "subject", "time",
        taxa = "B", verbose = FALSE
    )
    expect_identical(tr$results$taxon, "B")

    expect_error(
        mc_trends(tab, meta, "sample", "subject", "time",
            taxa = "C", verbose = FALSE
        ),
        "no row"
    )
})

test_that("a matrix input works and the taxon column is found either way", {
    meta <- mc_td()
    tab <- mc_tt(meta, list(X = function(d) 0.8 * d$time))
    m <- as.matrix(tab[, -1, drop = FALSE])
    rownames(m) <- tab$taxon

    from_df <- mc_trends(tab, meta, "sample", "subject", "time",
        verbose = FALSE
    )
    from_mat <- mc_trends(m, meta, "sample", "subject", "time",
        verbose = FALSE
    )
    expect_equal(from_df$results, from_mat$results)

    named <- tab
    names(named)[1] <- "OTU_ID"
    from_named <- mc_trends(named, meta, "sample", "subject", "time",
        taxon_col = "OTU_ID", verbose = FALSE
    )
    expect_equal(from_df$results, from_named$results)
})

test_that("out_dir writes the table and a log naming the model", {
    meta <- mc_td()
    tab <- mc_tt(meta, list(X = function(d) 0.8 * d$time))
    dir <- withr::local_tempdir()

    tr <- mc_trends(tab, meta, "sample", "subject", "time",
        out_dir = dir, verbose = FALSE
    )
    expect_setequal(basename(tr$files), "trends.tsv")
    back <- utils::read.delim(
        file.path(dir, "trends.tsv"),
        check.names = FALSE
    )
    expect_identical(back$taxon, tr$results$taxon)

    # The model that produced the table is on the object, since no log is
    # written for it to be read out of.
    expect_match(tr$settings$formula, "value ~ s(time, k =", fixed = TRUE)
    expect_identical(tr$settings$method, "REML")
    expect_match(tr$settings$test, "s(time)", fixed = TRUE)

    cm <- mc_compare(tab, meta, "sample", "subject", "time",
        group_col = "arm", out_dir = dir, verbose = FALSE
    )
    expect_setequal(basename(cm$files), "compare.tsv")
    expect_false(any(grepl("[.]txt$", list.files(dir))))
    expect_true(nzchar(cm$settings$reference))
    expect_output(print(cm), "test         :")
})

test_that("a trend run writes its table and no figures", {
    meta <- mc_td()
    tab <- mc_tt(meta, list(
        Diverges = function(d) ifelse(d$arm == "armB", 1.2 * d$time, 0)
    ))
    dir <- withr::local_tempdir()

    cm <- mc_compare(tab, meta, "sample", "subject", "time",
        group_col = "arm", out_dir = dir, verbose = FALSE
    )
    # Drawing is mc_trend_plot()'s job, so a run leaves nothing to look at
    # before its q-values have been read.
    expect_setequal(basename(cm$files), "compare.tsv")
    expect_length(list.files(dir, pattern = "[.]png$", recursive = TRUE), 0)

    imgs <- mc_trend_plot(cm, out_dir = dir, dpi = 72, verbose = FALSE)
    files <- list.files(dir, pattern = "[.]png$")
    # The name ends in a short hash of the full lineage, because two
    # species of one genus share every character a truncated name would
    # keep and the second figure would otherwise overwrite the first.
    expect_length(files, 1)
    expect_match(files, "^Diverges_[0-9a-f]{8}[.]png$")
    expect_s3_class(imgs, "ggplot")
})

test_that("a refit reproduces the fit the run reported", {
    meta <- mc_td()
    tab <- mc_tt(meta, list(X = function(d) 0.8 * d$time))

    tr <- mc_trends(tab, meta, "sample", "subject", "time", verbose = FALSE)
    again <- mc_trend_refit(tr, "X")

    # Models are not stored, so the figures and any follow-up refit. That is
    # only sound if the refit is the same fit: gam() with REML uses no RNG.
    st <- summary(again$fit)$s.table
    expect_equal(unname(st["s(time)", "edf"]), tr$results$edf[1],
        tolerance = 1e-10
    )
    expect_equal(
        unname(st["s(time)", ncol(st)]), tr$results$p_value[1],
        tolerance = 1e-10
    )
})

test_that("the result keeps the data but not the models", {
    meta <- mc_td()
    tab <- mc_tt(meta, list(X = function(d) 0.8 * d$time))

    tr <- mc_trends(tab, meta, "sample", "subject", "time", verbose = FALSE)
    expect_null(tr$models)
    expect_identical(nrow(tr$design), nrow(meta))
    expect_identical(dim(tr$abundance), c(1L, nrow(meta)))
})

test_that("two species of one genus get two files, not one", {
    meta <- mc_td()
    # Names that agree for far longer than any truncation would keep, and
    # differ only at the very end. This is the real case: a lineage runs
    # past 80 characters before the species is reached.
    long <- paste0(
        "k__Bacteria|p__Actinobacteria|c__Actinomycetia|",
        "o__Bifidobacteriales|f__Bifidobacteriaceae|g__Bifidobacterium|",
        "s__Bifidobacterium_"
    )
    nm <- paste0(long, c("breve", "longum"))

    tab <- mc_tt(meta, stats::setNames(
        list(function(d) 0.9 * d$time, function(d) -0.9 * d$time), nm
    ))

    dir <- withr::local_tempdir()
    tr <- mc_trends(tab, meta, "sample", "subject", "time", verbose = FALSE)
    mc_trend_plot(tr, out_dir = dir, dpi = 72, verbose = FALSE)

    expect_equal(sum(tr$results$q_value < 0.05, na.rm = TRUE), 2)
    expect_length(list.files(dir), 2)
})

test_that("a vector of taxa gives one figure each, whatever their q value", {
    meta <- mc_td()
    tab <- mc_tt(meta, list(
        Rising = function(d) 0.9 * d$time,
        Flat = function(d) 0
    ))

    tr <- mc_trends(tab, meta, "sample", "subject", "time", verbose = FALSE)
    expect_gt(tr$results$q_value[tr$results$taxon == "Flat"], 0.05)

    # Flat is not significant and is drawn anyway: naming it is already a
    # decision that it is worth looking at.
    ps <- mc_trend_plot(tr, c("Rising", "Flat"))
    expect_length(ps, 2)
    expect_named(ps, c("Rising", "Flat"))
    expect_s3_class(ps[["Flat"]], "ggplot")

    # One name gives the figure itself, not a list of one.
    one <- mc_trend_plot(tr, "Flat")
    expect_s3_class(one, "ggplot")

    dir <- withr::local_tempdir()
    mc_trend_plot(tr, "Flat", out_dir = dir, dpi = 72, verbose = FALSE)
    expect_match(list.files(dir), "^Flat_[0-9a-f]{8}[.]png$")
})

test_that("mc_trend_plot refuses what it cannot draw", {
    meta <- mc_td()
    tab <- mc_tt(meta, list(Flat = function(d) 0))
    tr <- mc_trends(tab, meta, "sample", "subject", "time", verbose = FALSE)

    expect_error(mc_trend_plot(tr, "Nope"), "Not in this result")
    # Nothing significant, so the default set is empty and the call says so
    # rather than returning nothing.
    expect_error(mc_trend_plot(tr), "nothing to draw")
    expect_error(mc_trend_plot(tr, 1), "character vector")
    expect_error(mc_trend_plot(tr$results, "Flat"), "mc_trends")
})


# k: one size, or a set of candidates chosen from per taxon.

test_that("k = NULL takes the largest basis the design supports", {
    meta <- mc_td(n_t = 6)
    tab <- mc_tt(meta, list(X = function(d) 0.8 * d$time))

    tr <- mc_trends(tab, meta, "sample", "subject", "time", verbose = FALSE)
    expect_identical(tr$settings$k, 5L)
    expect_identical(unique(tr$results$k), 5L)
    expect_match(tr$settings$formula, "k = 5", fixed = TRUE)
})

test_that("a single k is used as given and reported per taxon", {
    meta <- mc_td()
    tab <- mc_tt(meta, list(X = function(d) 0.8 * d$time))

    tr <- mc_trends(tab, meta, "sample", "subject", "time",
        k = 4, verbose = FALSE
    )
    expect_identical(unique(tr$results$k), 4L)
    expect_match(tr$settings$formula, "k = 4", fixed = TRUE)
})

test_that("several k are candidates, and each taxon keeps one of them", {
    meta <- mc_td()
    tab <- mc_tt(meta, list(
        Straight = function(d) 0.8 * d$time,
        Wavy = function(d) 2 * sin(d$time)
    ))

    tr <- mc_trends(tab, meta, "sample", "subject", "time",
        k = c(3, 4, 5), verbose = FALSE
    )
    expect_true(all(tr$results$k %in% c(3L, 4L, 5L)))
    expect_identical(tr$settings$k, c(3L, 4L, 5L))
    # No single k was fitted, so the recorded formula must not claim one.
    expect_match(tr$settings$formula, "chosen per taxon", fixed = TRUE)
})

test_that("the variance rule prefers the smallest basis that loses nothing", {
    meta <- mc_td()
    # A straight line needs no wiggle room, so the smallest candidate
    # reaches the largest one's deviance and wins.
    tab <- mc_tt(meta, list(Straight = function(d) 0.8 * d$time))

    tr <- mc_trends(tab, meta, "sample", "subject", "time",
        k = c(3, 4, 5), k_select = "variance", verbose = FALSE
    )
    expect_identical(tr$results$k, 3L)
})

test_that("a candidate k is capped by the time points a taxon has", {
    meta <- mc_td(n_t = 4)
    tab <- mc_tt(meta, list(X = function(d) 0.8 * d$time))

    tr <- mc_trends(tab, meta, "sample", "subject", "time",
        k = c(3, 4), verbose = FALSE
    )
    expect_true(all(tr$results$k <= 3L))
})

test_that("the figure is refitted at the k that taxon was tested with", {
    meta <- mc_td()
    tab <- mc_tt(meta, list(X = function(d) 2 * sin(d$time)))

    tr <- mc_trends(tab, meta, "sample", "subject", "time",
        k = c(3, 5), verbose = FALSE
    )
    rf <- mc_trend_refit(tr, "X")
    expect_identical(
        as.integer(rf$fit$smooth[[1]]$bs.dim), as.integer(tr$results$k[1])
    )
})

test_that("an impossible k is refused rather than silently lowered", {
    meta <- mc_td()
    tab <- mc_tt(meta, list(X = function(d) 0.8 * d$time))
    one <- function(...) {
        mc_trends(tab, meta, "sample", "subject", "time", ..., verbose = FALSE)
    }

    expect_error(one(k = 2), "at least 3")
    expect_error(one(k = c(3, 2)), "at least 3")
    expect_error(one(k = 99), "exceeds")
    expect_error(one(k = c(3, 99)), "exceeds")
    expect_error(one(k = "4"), "NULL, a single number, or a vector")
    expect_error(one(k = c(3, 5), k_var = 0), "greater than 0")
    expect_error(one(k = c(3, 5), k_var = 1.5), "greater than 0")
})

test_that("a grouped figure is one panel, not two", {
    meta <- mc_td()
    tab <- mc_tt(meta, list(X = function(d) 0.8 * d$time))

    cm <- mc_compare(tab, meta, "sample", "subject", "time",
        group_col = "arm", verbose = FALSE
    )
    p <- mc_trend_plot(cm, taxa = "X", verbose = FALSE)
    expect_s3_class(p, "ggplot")
    expect_false(inherits(p, "patchwork"))
})


test_that("leave-one-out picks a candidate and reports it", {
    meta <- mc_td()
    tab <- mc_tt(meta, list(
        Straight = function(d) 0.8 * d$time,
        Wavy = function(d) 2 * sin(d$time)
    ))

    tr <- mc_trends(tab, meta, "sample", "subject", "time",
        k = c(3, 4, 5), k_select = "loo", verbose = FALSE
    )
    expect_true(all(tr$results$k %in% c(3L, 4L, 5L)))
    expect_true(all(tr$results$converged))
})

test_that("the leave-one-out shortcut is the refitting answer", {
    # The closed form is only worth having if it agrees with the n refits
    # it replaces, so it is checked against them rather than against
    # itself. Smoothing parameters are held where REML put them, which is
    # the approximation the shortcut makes and the one being tested.
    skip_if_not_installed("mgcv")
    meta <- mc_td(n_sub = 8, n_t = 5)
    tab <- mc_tt(meta, list(X = function(d) 2 * sin(d$time)))
    d <- data.frame(
        value = as.numeric(unlist(tab[1, -1])),
        time = meta$time, subject = factor(meta$subject)
    )

    f <- value ~ s(time, k = 4) + s(subject, bs = "re")
    m <- mgcv::gam(f, data = d, method = "REML")
    shortcut <- mc_trend_loo_mse(m)

    brute <- mean(vapply(seq_len(nrow(d)), function(i) {
        mi <- mgcv::gam(f, data = d[-i, ], method = "REML", sp = m$sp)
        (d$value[i] - as.numeric(mgcv::predict.gam(mi, d[i, ])))^2
    }, numeric(1)))

    expect_equal(shortcut, brute, tolerance = 0.01)
})

test_that("an unknown k_select is refused", {
    meta <- mc_td()
    tab <- mc_tt(meta, list(X = function(d) 0.8 * d$time))
    expect_error(
        mc_trends(tab, meta, "sample", "subject", "time",
            k = c(3, 4), k_select = "bogus", verbose = FALSE
        ),
        "should be one of"
    )
})
