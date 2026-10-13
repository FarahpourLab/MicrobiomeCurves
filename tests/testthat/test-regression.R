# Numeric regression guard.
#
# These values pin the numbers the imputation returns for a fixed input,
# mask and seed.
#
# They were regenerated on 2026-10-12, when the estimator changed. Until then
# the package recomputed each subject's FPC scores with its own BLUP algebra;
# it now reads the scores fdapace itself used, from fp$xiEst, so the value is
# fdapace's own fitted curve at the nearest grid point. That makes the
# package agree bit for bit with the original analysis scripts
# (impute_FPCA_multi.R), which is what these numbers now represent. The
# previous constants came from the recomputation and no longer apply.
#
# These tests fail if a change to the FPCA, CI, clustering or outlier code
# alters the computed values. The published results depend on these numbers,
# so an intended change to the method requires regenerating the constants.
#
# Tolerance is 1e-6 rather than machine precision. FPCA uses a LAPACK
# eigendecomposition, and reference BLAS, OpenBLAS, MKL and Accelerate give
# slightly different results. Bioconductor builds on all three platforms, so a
# tighter tolerance would fail for reasons unrelated to this package. A real
# change to the method would exceed 1e-6; BLAS differences stay well below it.

test_that("imputed values for a fixed mask and seed have not changed", {
    expected_imputed <- c(
        -0.437956505400, -0.137865898166, 1.046200699461, -1.592710695464,
        1.387937360709, 0.811784415369, 0.849457710695, -1.835141107832,
        -1.048613426028, -0.506314640002, -0.550506957222, -0.930188319748
    )

    fit <- quiet_fit(demo_prep(), C = 1, use_outliers = TRUE, seed = 123)
    got <- mc_impute(fit)

    expect_equal(got$species, sprintf("Taxon%02d", 1:12))
    expect_equal(got$imputed_value, expected_imputed, tolerance = 1e-6)
})

test_that("the true values recovered at the masked cell are the input values", {
    fit <- quiet_fit(demo_prep(), C = 1, seed = 123)
    got <- mc_impute(fit)

    expect_equal(got$true_value, taxa_demo[["S01.3"]], tolerance = 1e-12)
})

test_that("squared error is consistent with true and imputed values", {
    fit <- quiet_fit(demo_prep(), C = 1, seed = 123)
    got <- mc_impute(fit)

    expect_equal(
        got$se, (got$true_value - got$imputed_value)^2, tolerance = 1e-12
    )
})

test_that("mc_metrics agrees with the long table it summarises", {
    fit <- quiet_fit(demo_prep(), C = 1, seed = 123)
    got <- mc_impute(fit)
    m <- mc_metrics(fit)

    expect_equal(m$MSE_overall, mean(got$se), tolerance = 1e-12)
    expect_equal(m$RMSE_overall, sqrt(m$MSE_overall), tolerance = 1e-12)
})

test_that("mc_run output for the bundled data has not changed", {
    run <- quiet_run(taxa_demo, taxon_col = "OTU_ID", C = 1, seed = 123)

    expect_equal(
        head(run$imputed$imputed_value, 6),
        c(0.692609755148, 0.565286316687, -0.688724509638,
            0.898977142979, 0.845289225360, -0.973382366093),
        tolerance = 1e-6
    )
    expect_equal(
        sum(run$imputed$imputed_value), 6.797706401138,
        tolerance = 1e-6
    )
})

test_that("results are reproducible across repeated runs", {
    a <- quiet_run(taxa_demo, taxon_col = "OTU_ID", C = 1, seed = 7)
    b <- quiet_run(taxa_demo, taxon_col = "OTU_ID", C = 1, seed = 7)

    expect_identical(a$imputed$imputed_value, b$imputed$imputed_value)
})

test_that("true values are NA when nothing was observed to compare against", {
    # mc_run() imputes data that was never observed, so there is no ground
    # truth and mc_metrics() does not give a meaningful score.
    run <- quiet_run(taxa_demo, taxon_col = "OTU_ID", C = 1)

    expect_true(all(is.na(run$imputed$true_value)))
    expect_true(all(is.na(run$imputed$se)))
})
