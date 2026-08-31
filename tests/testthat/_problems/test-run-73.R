# Extracted from test-run.R:73

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "MicrobiomeCurves", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
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
