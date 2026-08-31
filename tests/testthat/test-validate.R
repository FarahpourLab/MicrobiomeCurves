demo_bits <- function() {
    utils::data("taxa_demo", package = "MicrobiomeCurves", envir = environment())
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

test_that("a non-numeric sample column names itself", {
    b <- demo_bits()
    bad <- b$dat
    bad[[b$sub_cols[1]]] <- as.character(bad[[b$sub_cols[1]]])
    expect_error(
        mc_prepare(bad, "OTU_ID", b$mask),
        "not numeric.*S01\\.0"
    )
})

test_that("duplicated taxon names are rejected", {
    b <- demo_bits()
    dup <- rbind(b$dat, b$dat[1, ])
    expect_error(
        mc_prepare(dup, "OTU_ID", b$mask),
        "duplicated taxon names"
    )
})

test_that("a taxon with no data warns but does not stop the run", {
    b <- demo_bits()
    empty <- b$dat
    empty[1, b$sub_cols] <- NA_real_
    expect_warning(
        prep <- mc_prepare(empty, "OTU_ID", b$mask),
        "no observed values"
    )
    expect_s3_class(prep$mask_pairs, "data.frame")
})

test_that("a mask that matches nothing says what did not match", {
    b <- demo_bits()
    ghost <- data.frame(
        rep = "GHOST", time = b$times[1], stringsAsFactors = FALSE
    )
    expect_error(
        mc_prepare(b$dat, "OTU_ID", ghost),
        "subjects not in the table"
    )
    late <- data.frame(
        rep = b$subs[1], time = 999, stringsAsFactors = FALSE
    )
    expect_error(
        mc_prepare(b$dat, "OTU_ID", late),
        "time points not in the table"
    )
})

test_that("an empty mask is distinguished from a mismatched one", {
    b <- demo_bits()
    expect_error(
        mc_prepare(b$dat, "OTU_ID"),
        "marked no samples as missing"
    )
})

test_that("masking every sample is refused", {
    b <- demo_bits()
    all_masked <- matrix(0,
        nrow = length(b$subs), ncol = length(b$times),
        dimnames = list(b$subs, as.character(b$times))
    )
    expect_error(
        mc_prepare(b$dat, "OTU_ID", mask_matrix = all_masked),
        "every one of the .* samples as missing"
    )
})

test_that("a single subject is refused", {
    b <- demo_bits()
    one <- b$dat[, c("OTU_ID", grep(
        paste0("^", b$subs[1], "\\."), names(b$dat),
        value = TRUE
    ))]
    expect_error(
        mc_prepare(one, "OTU_ID", b$mask),
        "Only 1 subject was found"
    )
})

test_that("mask_list rejects elements that are not pairs", {
    b <- demo_bits()
    trap <- stats::setNames(list(b$times[2]), b$subs[1])
    expect_error(
        mc_prepare(b$dat, "OTU_ID", mask_list = trap),
        "must be a \\(subject, time\\) pair of length 2"
    )
    expect_error(
        mc_prepare(b$dat, "OTU_ID",
            mask_list = list(c(b$subs[1], "Tuesday"))
        ),
        "is not numeric"
    )
})

test_that("the documented mask_list forms both work", {
    b <- demo_bits()
    from_list <- mc_prepare(b$dat, "OTU_ID",
        mask_list = list(c(b$subs[1], b$times[2]))
    )
    from_df <- mc_prepare(b$dat, "OTU_ID", b$mask)
    expect_equal(from_list$mask_pairs, from_df$mask_pairs)
})

test_that("mc_metrics warns when there is no truth to score against", {
    b <- demo_bits()
    run <- suppressWarnings(
        run_wide_as_meta(b$dat, C = 1, verbose = FALSE)
    )
    expect_warning(
        m <- mc_metrics(run$fit),
        "No true values are available"
    )
    expect_true(is.nan(m$MSE_overall))
})

test_that("mc_detect_missing reports unusable columns", {
    b <- demo_bits()
    junk <- b$dat
    junk$free_text <- "note"
    expect_warning(
        mc_detect_missing(junk, "OTU_ID"),
        "do not match"
    )
})

test_that("a partly measured sample stops the run on every path", {
    b <- demo_bits()
    partial <- b$dat
    # One column observed for some taxa and NA for the rest: neither a
    # missing sample nor a complete one.
    partial[seq_len(3), b$sub_cols[1]] <- NA_real_

    expect_error(
        mc_prepare(partial, "OTU_ID", b$mask),
        "only partly measured"
    )
    expect_error(
        mc_detect_missing(partial, "OTU_ID"),
        "only partly measured"
    )
    expect_error(
        run_wide_as_meta(partial, C = 1, verbose = FALSE),
        "only partly measured"
    )
})

test_that("the partly measured error names the columns and counts", {
    b <- demo_bits()
    partial <- b$dat
    partial[seq_len(3), b$sub_cols[1]] <- NA_real_

    err <- tryCatch(
        mc_prepare(partial, "OTU_ID", b$mask),
        error = function(e) conditionMessage(e)
    )
    expect_match(err, b$sub_cols[1], fixed = TRUE)
    expect_match(err, "3 of 12 taxa NA", fixed = TRUE)
})

test_that("a fully NA sample is still fine, and is imputed", {
    b <- demo_bits()
    whole <- b$dat
    whole[[b$sub_cols[1]]] <- NA_real_
    # An entirely missing sample is the case the method exists for, so this
    # must not be caught by the partly-measured check.
    expect_no_error(mc_detect_missing(whole, "OTU_ID"))
})

test_that("kmeans_fd is a real alternative, not an alias", {
    b <- demo_bits()
    prep <- mc_prepare(b$dat, "OTU_ID", b$mask)
    f_fpca <- suppressWarnings(mc_fit(prep, C = 2))
    f_kfd <- suppressWarnings(
        mc_fit(prep, C = 2, cluster_method = "kmeans_fd")
    )

    expect_identical(f_fpca$cluster_method, "fpca")
    expect_identical(f_kfd$cluster_method, "kmeans_fd")

    sig <- function(f) {
        vapply(f$clusters, function(x) {
            paste(as.integer(x), collapse = "")
        }, character(1))
    }
    # The two routes group at least one taxon differently; if this ever
    # holds for every taxon the dispatch has silently stopped working.
    expect_false(identical(unname(sig(f_fpca)), unname(sig(f_kfd))))
})

test_that("an unknown cluster_method is rejected", {
    b <- demo_bits()
    prep <- mc_prepare(b$dat, "OTU_ID", b$mask)
    expect_error(
        mc_fit(prep, C = 2, cluster_method = "nope"),
        "should be one of"
    )
})
