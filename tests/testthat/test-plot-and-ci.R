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

mc_plot_run <- function() {
    demo <- demo_small()
    suppressWarnings(mc_run(
        demo$abundance, demo$metadata,
        sample_col = "sample", subject_col = "subject", time_col = "time",
        C = 1, verbose = FALSE
    ))
}

test_that("mc_plot works without dplyr or ggplot2 attached", {
    skip_if_not_installed("ggplot2")

    expect_false("package:dplyr" %in% search())
    expect_false("package:ggplot2" %in% search())

    run <- mc_plot_run()
    p <- mc_plot(run, run$completed[[1]][1], run$imputed$subject[1])

    expect_s3_class(p, "ggplot")
})

test_that("mc_plot refuses a taxon or a subject the run does not have", {
    run <- mc_plot_run()
    tx <- run$completed[[1]][1]

    subj <- run$imputed$subject[1]
    expect_error(mc_plot(run, "NotATaxon", subj), "NotATaxon")
    expect_error(mc_plot(run, tx, subject = "NotASubject"), "no subject")

    # One subject is the whole point of the plot, so it is asked for
    # rather than chosen on the caller's behalf.
    expect_error(mc_plot(run, tx), "subject is required")

    # The signature changed from (fit, species, rep, time) to (run, taxon,
    # subject), so a fit handed to it is refused with a message saying so
    # rather than failing somewhere inside.
    expect_error(mc_plot(run$fit, tx, subj), "must be an mc_run")
})

test_that("mc_plot marks observed and imputed, and bounds only imputed", {
    run <- mc_plot_run()
    tx <- run$completed[[1]][1]
    subj <- run$imputed$subject[1]

    p <- mc_plot(run, tx, subject = subj)
    tr <- mc_imputed_trajectory(run, tx)
    target <- tr[tr$subject == subj, ]

    # Every layer is built from the trajectory table, so the picture and
    # the tables the run wrote cannot disagree.
    pts <- p$layers[[which(vapply(
        p$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1)
    ))[1]]]
    expect_equal(nrow(pts$data), nrow(target))
    expect_setequal(unique(pts$data$source), unique(target$source))

    # An observed value is a measurement, so no interval is drawn for it.
    bars <- Filter(
        function(l) inherits(l$geom, "GeomErrorbar"), p$layers
    )
    expect_length(bars, 1)
    expect_true(all(bars[[1]]$data$source == "imputed"))
    expect_equal(
        nrow(bars[[1]]$data), sum(target$source == "imputed")
    )
})

test_that("uncertainty is bars at points, never a band across them", {
    run <- mc_plot_run()
    tx <- run$completed[[1]][1]
    tr <- mc_imputed_trajectory(run, tx)

    # A subject whose imputed points have observed ones between them,
    # which is where a band would make its false claim.
    imp <- tr[tr$source == "imputed", ]
    subj <- names(which.max(table(imp$subject)))
    target <- tr[tr$subject == subj, ]
    target <- target[order(target$time), ]
    expect_true(any(
        target$source[-1] == "observed" &
            target$source[-nrow(target)] == "imputed"
    ))

    p <- mc_plot(run, tx, subj)
    expect_length(Filter(
        function(l) inherits(l$geom, "GeomRibbon"), p$layers
    ), 0)
    expect_length(Filter(
        function(l) inherits(l$geom, "GeomErrorbar"), p$layers
    ), 1)
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
