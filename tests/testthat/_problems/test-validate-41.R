# Extracted from test-validate.R:41

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "TaxaTimeImpute", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
demo_bits <- function() {
    utils::data("taxa_demo", package = "TaxaTimeImpute", envir = environment())
    d <- taxa_demo
    sub_cols <- grep("\\.[0-9]+$", names(d), value = TRUE)
    subs <- unique(sub("\\..*$", "", sub_cols))
    tms <- sort(unique(as.numeric(sub("^.*\\.", "", sub_cols))))
    list(
        dat = d, sub_cols = sub_cols, subs = subs, times = tms,
        mask = data.frame(
            rep = subs[1], time = tms[2], stringsAsFactors = FALSE
        )
    )
}

# test -------------------------------------------------------------------------
b <- demo_bits()
empty <- b$dat
empty[1, b$sub_cols] <- NA_real_
expect_warning(
        prep <- tti_prepare(empty, "OTU_ID", b$mask),
        "no observed values"
    )
