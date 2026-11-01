# mc_imputed_trajectory().
#
# The function assembles what the run already holds, so what is tested is
# that it assembles it faithfully: the values are the completed table's, the
# observed/imputed split is the run's own record rather than a guess from
# which columns were created, and the intervals are the ones already written
# to disk.

mc_traj_run <- function() {
    demo <- demo_small()
    suppressWarnings(mc_run(
        demo$abundance, demo$metadata,
        sample_col = "sample", subject_col = "subject", time_col = "time",
        C = 1, verbose = FALSE
    ))
}

# Every gap in the bundled cohort is a sample that has no column at all, so
# the run creates a name for it and run$metadata$imputed is TRUE. The other
# kind of gap -- a column that exists and holds NA for every taxon, a sample
# that was planned and came back empty -- carries the caller's own name and
# is flagged FALSE. One is made here, because the distinction the test below
# is about only exists when both kinds are present.
mc_traj_run_both_kinds <- function() {
    demo <- demo_small()
    ab <- demo$abundance

    # Blank a sample belonging to a subject that can spare one and still
    # have two left, so the design does not refuse it.
    per <- table(demo$metadata$subject)
    donor <- names(per)[per >= 3][1]
    j <- which(demo$metadata$subject == donor)[1]
    ab[, j] <- NA_real_

    suppressWarnings(mc_run(
        ab, demo$metadata,
        sample_col = "sample", subject_col = "subject", time_col = "time",
        C = 1, verbose = FALSE
    ))
}

test_that("one row per subject-timepoint, with the promised columns", {
    run <- mc_traj_run()
    tx <- run$completed[[1]][1]
    tr <- mc_imputed_trajectory(run, tx)

    expect_identical(
        names(tr),
        c("taxon", "subject", "sample", "time", "time_label", "value",
          "source", "lower", "upper")
    )
    expect_identical(nrow(tr), nrow(run$metadata))
    expect_true(all(tr$taxon == tx))
    expect_setequal(unique(tr$source), c("observed", "imputed"))
    expect_false(anyNA(tr$value))
})

test_that("the values are the completed table's, cell for cell", {
    run <- mc_traj_run()
    tx <- run$completed[[1]][2]
    tr <- mc_imputed_trajectory(run, tx)

    i <- match(tx, run$completed[[1]])
    expect_equal(
        tr$value,
        as.numeric(run$completed[i, tr$sample]),
        tolerance = 0
    )
})

test_that("source comes from the run's record, not from created columns", {
    run <- mc_traj_run_both_kinds()
    tx <- run$completed[[1]][1]
    tr <- mc_imputed_trajectory(run, tx)

    imp <- run$imputed[run$imputed$species == tx, ]
    key <- function(s, t) paste(s, t, sep = "\r")
    got <- tr[tr$source == "imputed", ]
    expect_setequal(
        key(got$subject, got$time),
        key(as.character(imp$subject), imp$time)
    )

    # The distinction that matters: a column that existed but held NA for
    # every taxon was imputed, yet run$metadata$imputed is FALSE for it. A
    # flag read off the metadata would miss exactly those cells.
    created <- run$metadata$sample[run$metadata$imputed]
    imputed_samples <- tr$sample[tr$source == "imputed"]
    expect_true(length(setdiff(imputed_samples, c(created, NA))) > 0)
})

test_that("intervals match what the run already reported", {
    run <- mc_traj_run()
    tx <- run$completed[[1]][1]
    tr <- mc_imputed_trajectory(run, tx)
    unc <- suppressWarnings(mc_uncertainty(run, tx))

    imp <- tr[tr$source == "imputed", ]
    imp <- imp[order(imp$subject, imp$time), ]
    unc <- unc[order(unc$subject, unc$time), ]

    expect_equal(imp$value, unc$imputed, tolerance = 1e-12)
    expect_equal(imp$lower, unc$lower, tolerance = 1e-12)
    expect_equal(imp$upper, unc$upper, tolerance = 1e-12)
})

test_that("no se is offered, because the run does not store one", {
    run <- mc_traj_run()
    tr <- mc_imputed_trajectory(run, run$completed[[1]][1])
    expect_false("se" %in% names(tr))

    # The reason: run$imputed has the column but never fills it. Were that
    # to change, this function should start reporting it.
    expect_true(all(is.na(run$imputed$se)))
})

test_that("observed rows carry no interval", {
    run <- mc_traj_run()
    tr <- mc_imputed_trajectory(run, run$completed[[1]][1])
    obs <- tr[tr$source == "observed", ]

    expect_true(all(is.na(obs$lower)))
    expect_true(all(is.na(obs$upper)))
})

test_that("the interval brackets the imputed value", {
    run <- mc_traj_run()
    tr <- mc_imputed_trajectory(run, run$completed[[1]][1])
    imp <- tr[tr$source == "imputed" & !is.na(tr$lower), ]

    expect_gt(nrow(imp), 0)
    expect_true(all(imp$lower <= imp$value))
    expect_true(all(imp$upper >= imp$value))
    # The interval is symmetric, so the midpoint is the value itself.
    expect_equal((imp$lower + imp$upper) / 2, imp$value, tolerance = 1e-8)
})

test_that("every row names the column its value sits in", {
    run <- mc_traj_run()
    tr <- mc_imputed_trajectory(run, run$completed[[1]][1])

    # Including the columns the run had to create for a subject-timepoint
    # that had no sample of its own. That is what makes the value findable
    # in run$completed, and it implies why the cell was missing without a
    # column of its own saying so.
    expect_false(anyNA(tr$sample))
    expect_true(all(tr$sample %in% names(run$completed)))
    created <- run$metadata$sample[run$metadata$imputed]
    expect_true(all(created %in% tr$sample[tr$source == "imputed"]))
})

test_that("several taxa come back stacked and ordered", {
    run <- mc_traj_run()
    tx <- run$completed[[1]][1:3]
    tr <- mc_imputed_trajectory(run, tx)

    expect_identical(nrow(tr), 3L * nrow(run$metadata))
    expect_identical(unique(tr$taxon), tx)
    # Within a taxon, subject then time.
    one <- tr[tr$taxon == tx[2], ]
    expect_identical(one$time[order(
        match(one$subject, unique(one$subject)), one$time
    )], one$time)
})

test_that("each taxon's rows are the same as asking for it alone", {
    run <- mc_traj_run()
    tx <- run$completed[[1]][1:2]
    both <- mc_imputed_trajectory(run, tx)

    for (t1 in tx) {
        alone <- mc_imputed_trajectory(run, t1)
        rownames(alone) <- NULL
        part <- both[both$taxon == t1, ]
        rownames(part) <- NULL
        expect_equal(part, alone)
    }
})

test_that("a taxon that is not there is refused, and so is a repeat", {
    run <- mc_traj_run()

    expect_error(mc_imputed_trajectory(run, "NoSuchTaxon"), "No taxon called")
    expect_error(mc_imputed_trajectory(run, character(0)), "at least one")
    expect_error(
        mc_imputed_trajectory(run, rep(run$completed[[1]][1], 2)),
        "more than once"
    )
})

test_that("something that is not a run is refused", {
    msg <- "must be an object returned"
    expect_error(mc_imputed_trajectory(list(), "x"), msg)
    expect_error(mc_imputed_trajectory(mtcars, "x"), msg)
})

test_that("the empty-imputation guard leaves every row observed", {
    run <- mc_traj_run()
    tx <- run$completed[[1]][1]
    tr <- mc_imputed_trajectory(run, tx)
    base <- tr
    base$source <- "observed"
    base$lower <- NA_real_
    base$upper <- NA_real_

    # mc_run() refuses a table with nothing missing, so a run with no
    # imputed cells cannot be built through it. The guard is exercised
    # directly: handed an empty long table, nothing is marked.
    empty <- run$imputed[0, , drop = FALSE]
    stripped <- run
    stripped$imputed <- empty
    out <- MicrobiomeCurves:::mc_trajectory_mark(base, stripped, tx)

    expect_true(all(out$source == "observed"))
    expect_true(all(is.na(out$lower)))
    expect_true(all(is.na(out$upper)))

    stripped$imputed <- NULL
    out2 <- MicrobiomeCurves:::mc_trajectory_mark(base, stripped, tx)
    expect_true(all(out2$source == "observed"))
})
