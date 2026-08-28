# These tests run in a session where dplyr and patchwork are NOT attached,
# which is the condition that used to break mc_plot():
#
#   * `arrange()` was called bare, with no importFrom(dplyr, arrange)
#   * `p1 + p2` needs patchwork's method, but patchwork sat in Imports with no
#     NAMESPACE import, so its namespace never loaded
#
# Both failed for any user who had not run library(dplyr) / library(patchwork)
# first. These tests fail again if either import is dropped.

test_that("mc_plot works without dplyr or patchwork attached", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("patchwork")

    expect_false("package:dplyr" %in% search())
    expect_false("package:patchwork" %in% search())

    fit <- quiet_fit(demo_prep(), K = 1, seed = 123)

    p <- suppressWarnings(
        mc_plot(
            fit = fit,
            species_name = "Taxon01",
            rep_id = "S01",
            time_id = 3,
            ci_method = "analytic"
        )
    )

    expect_s3_class(p, "ggplot")
})

test_that("mc_plot rejects a species it never fitted", {
    fit <- quiet_fit(demo_prep(), K = 1, seed = 123)

    expect_error(
        suppressWarnings(
            mc_plot(
                fit,
                species_name = "NotATaxon", rep_id = "S01", time_id = 3
            )
        ),
        "clustering"
    )
})

test_that("analytic confidence intervals are well formed", {
    fit <- quiet_fit(demo_prep(), K = 1, seed = 123)

    ci <- suppressWarnings(
        mc_ci(
            fit = fit,
            species_name = "Taxon01",
            rep_id = "S01",
            time_id = 3,
            method = "analytic"
        )
    )

    expect_true("analytic" %in% names(ci))
    expect_true(all(c("time", "lower", "upper") %in% names(ci$analytic)))
    expect_gt(nrow(ci$analytic), 0)
    expect_true(all(ci$analytic$lower <= ci$analytic$upper))
    expect_false(anyNA(ci$analytic$lower))
    expect_false(anyNA(ci$analytic$upper))
})

test_that("the imputed value falls inside its own analytic interval", {
    fit <- quiet_fit(demo_prep(), K = 1, seed = 123)
    imp <- mc_impute(fit)

    ci <- suppressWarnings(
        mc_ci(fit, species_name = "Taxon01", rep_id = "S01",
            time_id = 3, method = "analytic")
    )

    at <- ci$analytic[which.min(abs(ci$analytic$time - 3)), ]
    value <- imp$imputed_value[imp$species == "Taxon01"]

    expect_gte(value, at$lower)
    expect_lte(value, at$upper)
})
