# The table methods take `...` so the generic can dispatch, but they do not
# forward it. Without a guard, anything extra is dropped in silence, and a
# call written against an older argument name runs and quietly does
# something else.

test_that("the old name for the cluster count is refused, not ignored", {
    d <- mc_demo_data()

    expect_error(
        mc_run(
            d$counts, d$metadata,
            sample_col = "sample", subject_col = "subject",
            time_col = "time", K = 1, verbose = FALSE
        ),
        "`K` is now `C`"
    )
})

test_that("a misspelled argument is refused", {
    d <- mc_demo_data()

    expect_error(
        mc_run(
            d$counts, d$metadata,
            sample_col = "sample", subject_col = "subject",
            time_col = "time", subject_cols = "subject", verbose = FALSE
        ),
        "Unused argument"
    )
})

test_that("the guard passes anything the method actually takes", {
    expect_silent(mc_check_dots())
    expect_error(mc_check_dots(K = 1), "`K` is now `C`")
    expect_error(mc_check_dots(nonsense = 1), "nonsense")
})

test_that("a run with the current argument names is unaffected", {
    d <- mc_demo_data()

    run <- suppressWarnings(mc_run(
        d$counts, d$metadata,
        sample_col = "sample", subject_col = "subject", time_col = "time",
        C = 1, verbose = FALSE
    ))

    expect_s3_class(run, "mc_run")
    expect_equal(run$fit$C, 1)
})
