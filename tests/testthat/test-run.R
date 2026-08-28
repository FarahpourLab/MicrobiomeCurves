test_that("mc_run fills every missing cell and reports what it did", {
    run <- quiet_run(taxa_demo, taxon_col = "OTU_ID", K = 1)

    expect_s3_class(run, "mc_run")
    expect_equal(nrow(run$missing), 3)
    expect_equal(run$n_failed, 0)
    expect_equal(nrow(run$imputed), nrow(taxa_demo) * 3)
    expect_false(anyNA(run$imputed$imputed_value))
})

test_that("observed values are never overwritten", {
    # Imputation must not modify measured values.
    run <- quiet_run(taxa_demo, taxon_col = "OTU_ID", K = 1)

    observed_cols <- setdiff(names(taxa_demo), c("OTU_ID", "S02.1", "S07.4"))

    for (cl in observed_cols) {
        expect_identical(
            run$completed[[cl]], taxa_demo[[cl]],
            info = paste("column", cl, "was modified")
        )
    }
    expect_identical(run$completed$taxon, taxa_demo$OTU_ID)
})

test_that("the completed table keeps the input layout", {
    run <- quiet_run(taxa_demo, taxon_col = "OTU_ID", K = 1)

    # Same rows; every observed sample comes back under its own name, and
    # the absent one is added under a name built from its subject and time.
    expect_equal(nrow(run$completed), nrow(taxa_demo))
    observed <- setdiff(names(taxa_demo), "OTU_ID")
    expect_true(all(observed %in% names(run$completed)))
    expect_true("S04_2" %in% names(run$completed))
    expect_false("S04_2" %in% names(taxa_demo))
})

test_that("previously missing cells come back populated", {
    run <- quiet_run(taxa_demo, taxon_col = "OTU_ID", K = 1)

    for (cl in c("S02.1", "S07.4", "S04_2")) {
        expect_false(anyNA(run$completed[[cl]]), info = paste("column", cl))
        expect_type(run$completed[[cl]], "double")
    }
})

test_that("mc_run refuses a table with nothing missing", {
    dat <- data.frame(OTU_ID = c("T1", "T2"))
    for (s in c("S1", "S2")) {
        for (tt in 0:3) dat[[paste0(s, ".", tt)]] <- c(1, 2)
    }

    expect_error(quiet_run(dat, taxon_col = "OTU_ID"), "nothing to impute")
})

test_that("mc_run warns about subjects with too few observations", {
    # min_observed = 3 makes every subject here "thin" without pushing any of
    # them below what the fitting code can handle.
    dat <- data.frame(OTU_ID = c("T1", "T2", "T3"))
    for (s in c("S1", "S2", "S3", "S4")) {
        for (tt in 0:3) dat[[paste0(s, ".", tt)]] <- rnorm(3)
    }
    dat[["S1.3"]] <- NA_real_

    w <- character()
    withCallingHandlers(
        run_wide_as_meta(dat, K = 1, min_observed = 4,
            verbose = FALSE),
        warning = function(x) {
            w <<- c(w, conditionMessage(x))
            invokeRestart("muffleWarning")
        }
    )

    expect_true(any(grepl("fewer than 4 observed time points", w)))
})

test_that("fitting survives a subject with too few observations", {
    # Regression guard for two crashes inherited from the original package:
    #
    #   use_outliers = TRUE -> detect_outliers_depth() returned one entry per
    #     subject FPCA could use, while mc_fit() named it against ALL subjects
    #     ("'names' attribute [4] must be the same length as the vector [3]").
    #
    #   either mode -> with no observations left for the target subject,
    #     mc_analytic_ci() built obs_idx via sapply() over an empty vector,
    #     giving an empty list, so `phi[obs_idx, , drop = FALSE]` failed
    #     ("invalid subscript type 'list'").
    #
    # The benchmark drivers never hit either, because they force subjects with
    # < 2 observations to fully observed. mc_run() does not, so users can.
    set.seed(1)
    dat <- data.frame(OTU_ID = c("T1", "T2", "T3"))
    for (s in c("S1", "S2", "S3", "S4")) {
        for (tt in 0:3) dat[[paste0(s, ".", tt)]] <- rnorm(3)
    }
    # S1 is left with a single observed point, which the mask then removes
    dat[["S1.0"]] <- NA_real_
    dat[["S1.1"]] <- NA_real_
    dat[["S1.2"]] <- NA_real_

    prep <- mc_prepare(
        dat, taxon_col = "OTU_ID",
        mask_list = data.frame(rep = c("S1", "S2"), time = c(3, 1))
    )

    for (uo in c(TRUE, FALSE)) {
        fit <- expect_no_error(
            suppressWarnings(mc_fit(prep, K = 1, use_outliers = uo))
        )
        got <- mc_impute(fit)
        expect_equal(nrow(got), nrow(dat) * 2)
        expect_true(all(is.finite(got$imputed_value)))
    }
})

test_that("a subject with no observations falls back to the population mean", {
    # With nothing subject-specific to condition on, the BLUP scores are zero
    # and the prediction is the fitted mean curve. Every subject that has lost
    # its whole trajectory must therefore receive the same value for a taxon.
    set.seed(2)
    dat <- data.frame(OTU_ID = c("T1", "T2"))
    for (s in c("S1", "S2", "S3", "S4", "S5")) {
        for (tt in 0:3) dat[[paste0(s, ".", tt)]] <- rnorm(2)
    }
    # S1 and S2 lose everything
    for (tt in 0:3) {
        dat[[paste0("S1.", tt)]] <- NA_real_
        dat[[paste0("S2.", tt)]] <- NA_real_
    }

    run <- quiet_run(dat, taxon_col = "OTU_ID", K = 1)
    imp <- run$imputed

    for (taxon in c("T1", "T2")) {
        for (tt in 0:3) {
            vals <- imp$imputed_value[
                imp$species == taxon & imp$time == tt &
                    imp$subject %in% c("S1", "S2")
            ]
            expect_equal(length(vals), 2)
            expect_equal(vals[1], vals[2], tolerance = 1e-10)
        }
    }
})

test_that("mc_run adds no logic of its own", {
    # The wrapper must return the same result as calling the three underlying
    # functions directly.
    run <- quiet_run(taxa_demo, taxon_col = "OTU_ID", K = 1, seed = 42)

    dat_full <- taxa_demo
    dat_full[["S04.2"]] <- NA_real_

    prep <- mc_prepare(
        dat = dat_full,
        taxon_col = "OTU_ID",
        mask_list = run$missing[, c("subject", "time")]
    )
    manual <- mc_impute(quiet_fit(prep, K = 1, seed = 42))

    expect_equal(run$imputed$imputed_value, manual$imputed_value)
    expect_equal(run$imputed$species, manual$species)
    expect_equal(run$imputed$subject, manual$rep)
    expect_equal(run$imputed$time, manual$time)
})

test_that("print.mc_run returns its input invisibly", {
    run <- quiet_run(taxa_demo, taxon_col = "OTU_ID", K = 1)

    expect_output(print(run), "MicrobiomeCurves run")
    expect_output(print(run), "cells imputed")
    expect_identical(withVisible(print(run))$visible, FALSE)
})
