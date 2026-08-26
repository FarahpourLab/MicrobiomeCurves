# Extracted from test-detect-missing.R:54

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "TaxaTimeImpute", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
dat <- data.frame(OTU_ID = c("T1", "T2", "T3"))
for (s in c("S1", "S2")) {
        for (tt in 0:3) dat[[paste0(s, ".", tt)]] <- c(1, 2, 3)
    }
dat[["S1.1"]][2] <- NA
info <- tti_detect_missing(dat, taxon_col = "OTU_ID")
