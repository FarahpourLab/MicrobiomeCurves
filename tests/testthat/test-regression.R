# Numeric regression guard.
#
# These values were produced by the original package (published as
# LongMiImpute) and verified to be bit-for-bit identical after the rename to
# TaxaTimeImpute. They pin the numbers the imputation returns for a fixed
# input, mask and seed.
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
        -0.438716933064, -0.139171897594, 1.055973932162, -1.608299264547,
        1.396902714826, 0.815203014745, 0.848228706263, -1.834692261179,
        -1.051379898447, -0.508832336939, -0.552737382916, -0.922789097269
    )

    fit <- quiet_fit(demo_prep(), K = 1, use_outliers = TRUE, seed = 123)
    got <- tti_impute(fit)

    expect_equal(got$species, sprintf("Taxon%02d", 1:12))
    expect_equal(got$imputed_value, expected_imputed, tolerance = 1e-6)
})

test_that("the true values recovered at the masked cell are the input values", {
    fit <- quiet_fit(demo_prep(), K = 1, seed = 123)
    got <- tti_impute(fit)

    expect_equal(got$true_value, taxa_demo[["S01.3"]], tolerance = 1e-12)
})

test_that("squared error is consistent with true and imputed values", {
    fit <- quiet_fit(demo_prep(), K = 1, seed = 123)
    got <- tti_impute(fit)

    expect_equal(
        got$se, (got$true_value - got$imputed_value)^2, tolerance = 1e-12
    )
})

test_that("tti_metrics agrees with the long table it summarises", {
    fit <- quiet_fit(demo_prep(), K = 1, seed = 123)
    got <- tti_impute(fit)
    m <- tti_metrics(fit)

    expect_equal(m$MSE_overall, mean(got$se), tolerance = 1e-12)
    expect_equal(m$RMSE_overall, sqrt(m$MSE_overall), tolerance = 1e-12)
})

test_that("tti_run output for the bundled data has not changed", {
    run <- quiet_run(taxa_demo, taxon_col = "OTU_ID", K = 1, seed = 123)

    expect_equal(
        head(run$imputed$imputed_value, 6),
        c(0.692524571439, 0.557688057361, -0.685291940706,
            0.899446132207, 0.837783251969, -0.977631605642),
        tolerance = 1e-6
    )
    expect_equal(
        sum(run$imputed$imputed_value), 6.773032458899,
        tolerance = 1e-6
    )
})

test_that("results are reproducible across repeated runs", {
    a <- quiet_run(taxa_demo, taxon_col = "OTU_ID", K = 1, seed = 7)
    b <- quiet_run(taxa_demo, taxon_col = "OTU_ID", K = 1, seed = 7)

    expect_identical(a$imputed$imputed_value, b$imputed$imputed_value)
})

test_that("true values are NA when nothing was observed to compare against", {
    # tti_run() imputes data that was never observed, so there is no ground
    # truth and tti_metrics() does not give a meaningful score.
    run <- quiet_run(taxa_demo, taxon_col = "OTU_ID", K = 1)

    expect_true(all(is.na(run$imputed$true_value)))
    expect_true(all(is.na(run$imputed$se)))
})
