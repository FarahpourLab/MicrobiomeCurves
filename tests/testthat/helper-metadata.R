# Drive the metadata interface from a table stored in the internal layout,
# so tests written around taxa_demo keep exercising the real entry point.
#
# The original column names are reused as sample identifiers, so a test can
# still refer to a column by the name it had in the fixture. Only the route
# into the package changes, not the data.
run_wide_as_meta <- function(dat, taxon_col = "OTU_ID", ...) {
    cols <- setdiff(names(dat), taxon_col)
    parsed <- TaxaTimeImpute:::tti_parse_cols(cols)

    mat <- as.matrix(dat[, cols, drop = FALSE])
    rownames(mat) <- dat[[taxon_col]]

    meta <- data.frame(
        sample = cols,
        subject = parsed$rep,
        time = parsed$time,
        stringsAsFactors = FALSE
    )

    tti_run(mat, meta,
        sample_col = "sample", subject_col = "subject", time_col = "time",
        ...
    )
}

