# Extracted from test-validate.R:159

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
partial <- b$dat
partial[seq_len(3), b$sub_cols[1]] <- NA_real_
expect_error(
        tti_prepare(partial, "OTU_ID", b$mask),
        "only partly measured"
    )
expect_error(
        tti_detect_missing(partial, "OTU_ID"),
        "only partly measured"
    )
expect_error(
        tti_run(partial, "OTU_ID", K = 1, verbose = FALSE),
        "only partly measured"
    )
