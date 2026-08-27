# Extracted from test-run.R:141

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "TaxaTimeImpute", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
set.seed(2)
dat <- data.frame(OTU_ID = c("T1", "T2"))
for (s in c("S1", "S2", "S3", "S4", "S5")) {
        for (tt in 0:3) dat[[paste0(s, ".", tt)]] <- rnorm(2)
    }
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
                    imp$rep %in% c("S1", "S2")
            ]
            expect_equal(length(vals), 2)
            expect_equal(vals[1], vals[2], tolerance = 1e-10)
        }
    }
