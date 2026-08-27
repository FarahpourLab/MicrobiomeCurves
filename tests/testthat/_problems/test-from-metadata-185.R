# Extracted from test-from-metadata.R:185

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "TaxaTimeImpute", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
study <- function(drop = NULL, blank = NULL) {
    subjects <- paste0("M", sprintf("%02d", 1:5))
    days <- c(0, 7, 14)

    meta <- expand.grid(
        subject = subjects, day = days,
        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
    meta$sample <- paste0(meta$subject, "_d", meta$day)
    meta <- meta[order(meta$subject, meta$day), c("sample", "subject", "day")]
    if (!is.null(drop)) meta <- meta[meta$sample != drop, ]
    rownames(meta) <- NULL

    taxa <- paste0("Taxon", sprintf("%02d", 1:4))
    set.seed(42)
    m <- matrix(
        rnorm(length(taxa) * nrow(meta)),
        nrow = length(taxa), dimnames = list(taxa, meta$sample)
    )
    if (!is.null(blank)) m[, blank] <- NA_real_

    list(
        abundance = data.frame(taxon = taxa, m, check.names = FALSE),
        metadata = meta, taxa = taxa, subjects = subjects, days = days
    )
}
build <- function(s, ...) {
    tti_from_metadata(
        s$abundance, s$metadata,
        sample_col = "sample", subject_col = "subject", time_col = "day",
        taxon_col = "taxon", verbose = FALSE, ...
    )
}

# test -------------------------------------------------------------------------
s <- study(drop = "M03_d7")
run <- suppressWarnings(tti_run(
        s$abundance,
        metadata = s$metadata, sample_col = "sample",
        subject_col = "subject", time_col = "day",
        taxon_col = "taxon", K = 1, verbose = FALSE
    ))
cols <- setdiff(names(run$completed), "taxon")
expect_equal(
        cols[10:12],
        c("M03_d0", "M03_7", "M03_d14")
    )
