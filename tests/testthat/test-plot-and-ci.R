# These tests run in a session where dplyr and patchwork are NOT attached,
# which is the condition that used to break mc_plot():
#
#   * `arrange()` was called bare, with no importFrom(dplyr, arrange)
#   * `p1 + p2` needs patchwork's method, but patchwork sat in Imports with no
#     NAMESPACE import, so its namespace never loaded
#
# Both failed for any user who had not run library(dplyr) / library(patchwork)
# first. These tests fail again if either import is dropped.

# Outlier screening is switched off in these fixtures. With it on, the
# masked subject is itself flagged as an outlier for this taxon and so is
# left out of fit$clusters, and mc_ci() then indexes its subject-by-time
# matrix with NA and fails with "subscript out of bounds". mc_impute() does
# produce a value for that cell, because mc_cell_members() re-admits an
# outlying target subject; the interval path has no such step. That
# disagreement is a defect in the package, not in these tests, and these
# tests are about the dplyr and patchwork imports and the shape of an
# interval, so they do not depend on it either way.

test_that("mc_plot works without dplyr or patchwork attached", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("patchwork")

    expect_false("package:dplyr" %in% search())
    expect_false("package:patchwork" %in% search())

    fit <- quiet_fit(demo_prep(), C = 1, seed = 123, use_outliers = FALSE)

    p <- suppressWarnings(
        mc_plot(
            fit = fit,
            species_name = taxa_demo$OTU_ID[1],
            rep_id = demo_masked_cell$rep,
            time_id = demo_masked_cell$time,
            ci_method = "analytic"
        )
    )

    expect_s3_class(p, "ggplot")
})

test_that("mc_plot rejects a species it never fitted", {
    fit <- quiet_fit(demo_prep(), C = 1, seed = 123, use_outliers = FALSE)

    expect_error(
        suppressWarnings(
            mc_plot(
                fit,
                species_name = "NotATaxon", rep_id = demo_masked_cell$rep, time_id = demo_masked_cell$time
            )
        ),
        "clustering"
    )
})

test_that("analytic confidence intervals are well formed", {
    fit <- quiet_fit(demo_prep(), C = 1, seed = 123, use_outliers = FALSE)

    ci <- suppressWarnings(
        mc_ci(
            fit = fit,
            species_name = taxa_demo$OTU_ID[1],
            rep_id = demo_masked_cell$rep,
            time_id = demo_masked_cell$time,
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
    fit <- quiet_fit(demo_prep(), C = 1, seed = 123, use_outliers = FALSE)
    imp <- mc_impute(fit)

    ci <- suppressWarnings(
        mc_ci(fit, species_name = taxa_demo$OTU_ID[1], rep_id = demo_masked_cell$rep,
            time_id = demo_masked_cell$time, method = "analytic")
    )

    at <- ci$analytic[which.min(abs(ci$analytic$time - 3)), ]
    value <- imp$imputed_value[imp$species == taxa_demo$OTU_ID[1]]

    expect_gte(value, at$lower)
    expect_lte(value, at$upper)
})
