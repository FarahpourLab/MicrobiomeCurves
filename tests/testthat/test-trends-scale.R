# normalization and transform in mc_trends() and mc_compare().
#
# The vocabulary is MaAsLin 2's, so these tests pin the two things that
# matter: that each option does what that package's own code does, and that
# the one option meant to agree with this package's imputation path agrees
# with it exactly rather than nearly. The axes are also tested, because TSS
# works down a sample and LOG along a taxon, and swapping them would be a
# silent wrong answer rather than an error.

mc_scale_meta <- function(n_sub = 12, n_t = 6) {
    d <- expand.grid(
        time = seq_len(n_t), subject = paste0("S", seq_len(n_sub)),
        stringsAsFactors = FALSE
    )
    d <- d[order(d$subject, d$time), ]
    d$sample <- paste0(d$subject, "_T", d$time)
    half <- ceiling(n_sub / 2)
    d$arm <- ifelse(
        as.integer(sub("^S", "", d$subject)) <= half, "armA", "armB"
    )
    rownames(d) <- NULL
    d
}

mc_abund <- function(n_sub = 12, n_t = 6, n_tax = 4, seed = 7) {
    meta <- mc_scale_meta(n_sub = n_sub, n_t = n_t)
    m <- withr::with_seed(seed, {
        matrix(
            stats::rpois(n_tax * nrow(meta), lambda = 40),
            nrow = n_tax,
            dimnames = list(paste0("Taxon", seq_len(n_tax)), meta$sample)
        )
    })
    # Zeros in one taxon only, so the two zero rules can be told apart.
    m["Taxon2", seq_len(5)] <- 0
    list(meta = meta, mat = m,
        tab = data.frame(taxon = rownames(m), m, check.names = FALSE))
}


test_that("the default leaves the table exactly as supplied", {
    d <- mc_abund()
    tr <- mc_trends(d$tab, d$meta, "sample", "subject", "time",
        verbose = FALSE
    )

    expect_identical(tr$settings$normalization, "NONE")
    expect_identical(tr$settings$transform, "NONE")
    # The fitted matrix is the input, not a rescaling of it. This is what
    # keeps mc_trends(run$completed, ...) correct, since that table is
    # already CLR and any normalisation would be a second transform.
    expect_equal(
        tr$abundance,
        d$mat[, match(colnames(tr$abundance), colnames(d$mat))]
    )
})


test_that("CLR in mc_trends is the same transform the imputation uses", {
    d <- mc_abund()

    # The property the argument exists for: a trend test on an abundance
    # and an imputation of that same table stand on identical values, so
    # comparing the two measures the imputed cells and nothing else.
    tr <- mc_trends(d$tab, d$meta, "sample", "subject", "time",
        normalization = "CLR", verbose = FALSE
    )
    direct <- mc_clr(d$mat, "auto")

    expect_equal(
        tr$abundance,
        direct[, match(colnames(tr$abundance), colnames(direct))]
    )
    expect_identical(tr$settings$normalization, "CLR")
    expect_identical(tr$settings$pseudocount, "auto")
})


test_that("a numeric pseudocount is added everywhere, not just to zeros", {
    d <- mc_abund()

    auto <- mc_trends(d$tab, d$meta, "sample", "subject", "time",
        normalization = "CLR", verbose = FALSE
    )
    one <- mc_trends(d$tab, d$meta, "sample", "subject", "time",
        normalization = "CLR", pseudocount = 1, verbose = FALSE
    )

    expect_equal(one$abundance, mc_clr(d$mat, 1)[, colnames(one$abundance)])
    expect_false(isTRUE(all.equal(auto$abundance, one$abundance)))
    expect_identical(one$settings$pseudocount, 1)
})


test_that("pseudocount = \"auto\" uses the package's shared multiplier", {
    x <- c(0, 4, 6, 10)
    expect_equal(
        mc_replace_zeros(x, "auto"),
        c(4 * mc_zero_mult, 4, 6, 10)
    )
    # Named rather than written out, because mc_raw_scale_info() has to undo
    # exactly this replacement to decide which imputed values sit below the
    # detection limit.
    expect_equal(mc_zero_mult, 0.65)
})


test_that("an unusable pseudocount is refused before any fitting", {
    d <- mc_abund()
    expect_error(
        mc_trends(d$tab, d$meta, "sample", "subject", "time",
            normalization = "CLR", pseudocount = -1, verbose = FALSE
        ),
        "positive number"
    )
})


test_that("CLR cannot be combined with any transform", {
    d <- mc_abund()
    for (tf in c("LOG", "LOGIT", "AST")) {
        expect_error(
            mc_trends(d$tab, d$meta, "sample", "subject", "time",
                normalization = "CLR", transform = tf, verbose = FALSE
            ),
            "cannot be combined"
        )
    }
    expect_silent(mc_trend_check_scale("CLR", "NONE"))
})


test_that("TSS works down a sample and makes each one sum to one", {
    d <- mc_abund()
    tss <- mc_trend_tss(d$mat)

    expect_equal(colSums(tss), rep(1, ncol(d$mat)),
        ignore_attr = TRUE
    )
    # Per sample, so a taxon's values are not comparable across samples
    # before TSS and are after it. The row sums carry no such constraint.
    expect_false(isTRUE(all.equal(
        unname(rowSums(tss)), rep(1, nrow(d$mat))
    )))
})


test_that("TSS leaves an empty sample as NA rather than NaN", {
    m <- matrix(c(1, 2, 0, 0), nrow = 2,
        dimnames = list(c("a", "b"), c("s1", "s2"))
    )
    tss <- mc_trend_tss(m)

    expect_equal(tss[, "s1"], c(a = 1 / 3, b = 2 / 3))
    # NaN would survive every is.na() guard downstream and reach the GAM;
    # NA is dropped from the fit the same way a gap is.
    expect_true(all(is.na(tss[, "s2"])))
    expect_false(any(is.nan(tss[, "s2"])))
})


test_that("LOG works along a taxon, with MaAsLin 2's halved minimum", {
    m <- matrix(c(0, 2, 4, 100, 200, 400), nrow = 2, byrow = TRUE,
        dimnames = list(c("sparse", "dense"), c("s1", "s2", "s3"))
    )
    lg <- mc_trend_log(m)

    # The replacement is half of the smallest positive value in that taxon,
    # which is 2 for "sparse" and not 100, the study-wide minimum. MaAsLin 2
    # does apply(features, 2, LOG) on samples-by-features, so its axis is
    # the taxon too.
    expect_equal(lg["sparse", ], c(s1 = log2(1), s2 = log2(2), s3 = log2(4)))
    expect_equal(lg["dense", ], log2(c(s1 = 100, s2 = 200, s3 = 400)))
})


test_that("a taxon that is zero everywhere leaves LOG as NA, not an error", {
    m <- matrix(c(0, 0, 0, 1, 2, 3), nrow = 2, byrow = TRUE,
        dimnames = list(c("absent", "present"), c("s1", "s2", "s3"))
    )
    lg <- mc_trend_log(m)

    expect_true(all(is.na(lg["absent", ])))
    expect_equal(lg["present", ], log2(c(s1 = 1, s2 = 2, s3 = 3)))
})


test_that("AST matches MaAsLin 2 and refuses values off its domain", {
    m <- matrix(c(0, 0.25, 1), nrow = 1,
        dimnames = list("a", c("s1", "s2", "s3"))
    )
    expect_equal(
        mc_trend_ast(m)["a", ],
        sign(c(0, 0.25, 1)) * asin(sqrt(abs(c(0, 0.25, 1)))),
        ignore_attr = TRUE
    )
    expect_error(mc_trend_ast(m * 2), "between -1 and 1")
})


test_that("LOGIT reports the boundary values it moves to the middle", {
    m <- matrix(c(0, 0.5, 1), nrow = 1,
        dimnames = list("a", c("s1", "s2", "s3"))
    )
    expect_warning(out <- mc_trend_logit(m), "exactly 0 or 1")

    # MaAsLin 2's rule, kept for comparability: an infinite logit becomes 0,
    # which is a proportion of one half. The warning exists because that is
    # the middle of the scale, not the end of it.
    expect_equal(out["a", ], c(s1 = 0, s2 = 0, s3 = 0))
    expect_error(mc_trend_logit(m * 3), "between 0 and 1")
})


test_that("the scale is applied to the whole table, not to taxa selected", {
    d <- mc_abund()

    few <- mc_trends(d$tab, d$meta, "sample", "subject", "time",
        taxa = c("Taxon1", "Taxon2"), normalization = "CLR", verbose = FALSE
    )
    all_tax <- mc_trends(d$tab, d$meta, "sample", "subject", "time",
        normalization = "CLR", verbose = FALSE
    )

    # CLR centres on the geometric mean over the whole composition, so
    # selecting two taxa first would change every value. The selection has
    # to happen after the transform.
    expect_equal(
        few$abundance,
        all_tax$abundance[c("Taxon1", "Taxon2"), colnames(few$abundance)]
    )
    expect_equal(
        few$results[order(few$results$taxon), "p_value"],
        all_tax$results[
            match(c("Taxon1", "Taxon2"), all_tax$results$taxon), "p_value"
        ]
    )
})


test_that("the scale reaches the settings and the printout", {
    d <- mc_abund()
    out <- withr::local_tempdir()

    tr <- mc_trends(d$tab, d$meta, "sample", "subject", "time",
        normalization = "TSS", transform = "LOG",
        out_dir = out, verbose = FALSE
    )

    expect_identical(tr$settings$normalization, "TSS")
    expect_identical(tr$settings$transform, "LOG")

    # A trend table with no record of the scale it was fitted on cannot be
    # reproduced, so the printout has to carry it. No log is written.
    expect_setequal(basename(tr$files), "trends.tsv")
    expect_false(any(grepl("[.]txt$", list.files(out))))
    expect_output(print(tr), "input scale  : TSS \\+ LOG")

    plain <- mc_trends(d$tab, d$meta, "sample", "subject", "time",
        verbose = FALSE
    )
    expect_false(any(grepl("input scale", capture.output(print(plain)))))
})


test_that("the CLR scale line names the multiplier actually used", {
    d <- mc_abund()

    msg <- capture.output(
        tr <- mc_trends(d$tab, d$meta, "sample", "subject", "time",
            normalization = "CLR"
        ),
        type = "message"
    )

    # The multiplier actually used is named where the run is reported, so
    # someone reading the session can see which rule the zeros took.
    expect_true(any(grepl("pseudocount \"auto\"", msg, fixed = TRUE)))
    expect_true(any(grepl(as.character(mc_zero_mult), msg, fixed = TRUE)))
    expect_identical(tr$settings$pseudocount, "auto")

    # One formatter behind the lines, whoever reads them.
    lines <- mc_trend_scale_lines("CLR", "NONE", "auto")
    expect_true(all(vapply(
        lines, function(l) any(grepl(l, msg, fixed = TRUE)), logical(1)
    )))
})


test_that("mc_compare takes the same two arguments as mc_trends", {
    d <- mc_abund()

    cm <- mc_compare(d$tab, d$meta, "sample", "subject", "time",
        group_col = "arm", normalization = "CLR", verbose = FALSE
    )

    expect_identical(cm$settings$normalization, "CLR")
    expect_identical(cm$settings$transform, "NONE")
    expect_equal(
        cm$abundance,
        mc_clr(d$mat, "auto")[, colnames(cm$abundance)]
    )
    expect_error(
        mc_compare(d$tab, d$meta, "sample", "subject", "time",
            group_col = "arm", normalization = "CLR", transform = "LOG",
            verbose = FALSE
        ),
        "cannot be combined"
    )
})


test_that("an unknown option is refused rather than ignored", {
    d <- mc_abund()
    expect_error(
        mc_trends(d$tab, d$meta, "sample", "subject", "time",
            normalization = "CSS", verbose = FALSE
        ),
        "should be one of"
    )
    expect_error(
        mc_trends(d$tab, d$meta, "sample", "subject", "time",
            transform = "SQRT", verbose = FALSE
        ),
        "should be one of"
    )
})


test_that("a CLR table refused by TSS or LOG says what to pass instead", {
    meta <- mc_scale_meta()
    tab <- withr::with_seed(1, {
        m <- matrix(stats::rnorm(nrow(meta), sd = 0.3), nrow = 1,
            dimnames = list("A", meta$sample)
        )
        data.frame(taxon = "A", m, check.names = FALSE)
    })

    expect_error(
        mc_trends(tab, meta, "sample", "subject", "time",
            normalization = "TSS", verbose = FALSE
        ),
        "non-negative"
    )
    expect_error(
        mc_trends(tab, meta, "sample", "subject", "time",
            transform = "LOG", verbose = FALSE
        ),
        "non-negative"
    )
})
