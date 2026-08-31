# Input validation for mc_prepare().
#
# Before these checks existed a malformed table reached the mask filter and
# died with "No valid masked pairs.", which named the wrong culprit: a
# character column, a duplicated taxon name and a mask whose subjects did not
# exist in the table all produced that same message. Each cause now reports
# itself.

#' Abbreviate a character vector for use in an error message
#'
#' @param x Character vector to render.
#' @param n Maximum number of elements to show before eliding.
#'
#' @return A single string of comma-separated values, with a trailing
#'   `"... (N more)"` when `x` is longer than `n`.
#'
#' @keywords internal
#' @noRd
mc_fmt_some <- function(x, n = 5L) {
    x <- as.character(x)
    if (length(x) <= n) {
        return(paste(x, collapse = ", "))
    }
    paste0(
        paste(x[seq_len(n)], collapse = ", "),
        "... (", length(x) - n, " more)"
    )
}

#' Check that every parsed sample column is numeric
#'
#' @param dat The user's data frame.
#' @param col_map Data frame of parsed columns, with a `col` column.
#'
#' @return `NULL`, invisibly. Called for the error it raises.
#'
#' @keywords internal
#' @noRd
mc_check_numeric_cols <- function(dat, col_map) {
    bad <- col_map$col[!vapply(dat[col_map$col], is.numeric, logical(1))]
    if (length(bad) > 0) {
        stop(
            "These sample columns are not numeric: ", mc_fmt_some(bad),
            ". Abundances must be numeric; a column read as text usually ",
            "means the file has a non-numeric placeholder for missing ",
            "values. Read it with na.strings= so those become NA.",
            call. = FALSE
        )
    }
    invisible(NULL)
}

#' Check that taxon identifiers are usable
#'
#' @param dat The user's data frame.
#' @param taxon_col Character name of the taxon column.
#'
#' @return `NULL`, invisibly. Errors on duplicates or missing identifiers.
#'
#' @keywords internal
#' @noRd
mc_check_taxa <- function(dat, taxon_col) {
    taxa <- dat[[taxon_col]]

    if (anyNA(taxa) || any(!nzchar(trimws(as.character(taxa))))) {
        stop(
            "Column '", taxon_col, "' has blank or NA taxon names. ",
            "Every row needs an identifier.",
            call. = FALSE
        )
    }

    dup <- unique(taxa[duplicated(taxa)])
    if (length(dup) > 0) {
        stop(
            "Column '", taxon_col, "' has duplicated taxon names: ",
            mc_fmt_some(dup),
            ". Results are reported per taxon name, so names must be ",
            "unique. Aggregate the duplicates or make the names distinct.",
            call. = FALSE
        )
    }
    invisible(NULL)
}

#' Check that the design can support a functional fit
#'
#' @description
#' FPCA borrows strength across subjects, so a single subject can never be
#' fitted. This is reported here rather than surfacing later as an
#' unexplained column of `NA` values.
#'
#' @param reps Character vector of subject identifiers found in the table.
#'
#' @return `NULL`, invisibly. Errors when fewer than two subjects are present.
#'
#' @keywords internal
#' @noRd
mc_check_design <- function(reps) {
    if (length(reps) >= 2) {
        return(invisible(NULL))
    }

    named <- if (length(reps) == 1) c(" ('", reps, "')") else ""
    stop(
        "Only ", length(reps), " subject was found", named,
        ". Imputation borrows information across subjects, so at least ",
        "two are required. Check that the sample columns are named ",
        "'<subject>.<time>'.",
        call. = FALSE
    )
}

#' Warn about taxa that carry no data
#'
#' @param dat The user's data frame.
#' @param taxon_col Character name of the taxon column.
#' @param col_map Data frame of parsed columns, with a `col` column.
#'
#' @return `NULL`, invisibly. Called for the warning it may raise.
#'
#' @keywords internal
#' @noRd
mc_warn_empty_taxa <- function(dat, taxon_col, col_map) {
    vals <- as.matrix(dat[col_map$col])
    empty <- rowSums(!is.na(vals)) == 0
    if (any(empty)) {
        warning(
            sum(empty), " taxa have no observed values at all: ",
            mc_fmt_some(dat[[taxon_col]][empty]),
            ". They cannot be imputed and will be returned as NA.",
            call. = FALSE
        )
    }
    invisible(NULL)
}

#' Explain why no masked pair survived the filter
#'
#' @description
#' Distinguishes a mask that names subjects or times absent from the table
#' from a mask that was empty to begin with, so the message points at the
#' actual mismatch.
#'
#' @param mask_pairs Data frame of requested pairs, before filtering.
#' @param reps Character vector of subject identifiers in the table.
#' @param times Numeric vector of time points in the table.
#'
#' @return This function always throws.
#'
#' @keywords internal
#' @noRd
mc_stop_no_mask_pairs <- function(mask_pairs, reps, times) {
    if (nrow(mask_pairs) == 0) {
        stop(
            "The mask marked no samples as missing. Supply at least one ",
            "(subject, time) pair, or leave the mask arguments NULL and use ",
            "mc_run() to impute the samples that are absent from the table.",
            call. = FALSE
        )
    }

    bad_rep <- setdiff(unique(mask_pairs$rep), reps)
    bad_time <- setdiff(unique(mask_pairs$time), times)

    detail <- character(0)
    if (length(bad_rep) > 0) {
        detail <- c(
            detail,
            "subjects not in the table: ", mc_fmt_some(bad_rep),
            " (the table has ", mc_fmt_some(reps), ")"
        )
    }
    if (length(bad_time) > 0) {
        if (length(detail) > 0) detail <- c(detail, "; ")
        detail <- c(
            detail,
            "time points not in the table: ", mc_fmt_some(bad_time),
            " (the table has ", mc_fmt_some(times), ")"
        )
    }

    if (length(detail) == 0) {
        detail <- c(
            "The requested (subject, time) pairs do not correspond to ",
            "any column."
        )
    } else {
        detail <- c("Mismatch: ", detail, ".")
    }

    stop("The mask matched no sample in the table. ", detail, call. = FALSE)
}

#' Run every check that applies to the table itself
#'
#' @param dat The user's data frame.
#' @param taxon_col Character name of the taxon column.
#' @param col_map Data frame of parsed columns.
#' @param reps Character vector of subject identifiers.
#'
#' @return `NULL`, invisibly. Called for its errors and warnings.
#'
#' @keywords internal
#' @noRd
mc_check_inputs <- function(dat, taxon_col, col_map, reps) {
    mc_check_numeric_cols(dat, col_map)
    mc_check_taxa(dat, taxon_col)
    mc_check_design(reps)
    mc_check_partial_na(dat, col_map)
    mc_warn_empty_taxa(dat, taxon_col, col_map)
    invisible(NULL)
}

#' Gather the requested mask into one table of pairs
#'
#' @param mask_list Optional list or data frame of (subject, time) pairs.
#' @param mask_matrix Optional matrix of observed/missing flags.
#'
#' @return Data frame with `rep` and `time`, possibly with no rows.
#'
#' @keywords internal
#' @noRd
mc_collect_mask_pairs <- function(mask_list, mask_matrix) {
    pairs <- tibble(rep = character(), time = numeric())

    if (!is.null(mask_list)) {
        pairs <- dplyr::bind_rows(
            pairs, mc_mask_pairs_from_list(mask_list)
        )
    }
    if (!is.null(mask_matrix)) {
        pairs <- dplyr::bind_rows(
            pairs, mc_mask_pairs_from_matrix(mask_matrix)
        )
    }
    pairs
}

#' Warn about columns that were skipped or partly empty
#'
#' @description
#' Shared by [mc_detect_missing()] and the [mc_run()] input check, so that
#' calling `mc_detect_missing()` directly reports these too. Both are also
#' returned in the result, but only as fields the caller has to think to look
#' at.
#'
#' @param unparsed Character vector of column names that did not parse.
#' @param partial_na Data frame of columns holding some but not only `NA`.
#'
#' @return `NULL`, invisibly. Called for its warnings.
#'
#' @keywords internal
#' @noRd
mc_warn_column_issues <- function(unparsed, partial_na) {
    if (length(unparsed) > 0) {
        warning(
            "Ignoring ", length(unparsed),
            " column(s) that do not match ",
            "'<subject>.<time>': ", mc_fmt_some(unparsed),
            call. = FALSE
        )
    }

    mc_stop_partial_na(partial_na$col, partial_na$n_na, partial_na$n_taxa)
    invisible(NULL)
}

#' Refuse a sample that is only partly measured
#'
#' @description
#' A sample is one column, holding a value for every taxon. The method imputes
#' whole missing samples, so a column is expected to be either fully observed
#' or fully `NA`. A column with some taxa present and others `NA` is neither,
#' and the scattered `NA` values would silently drop out of the fit. That is
#' reported as an error rather than a warning so the data can be corrected
#' before any result is produced.
#'
#' @param cols Character vector of offending column names.
#' @param n_na Integer vector, how many taxa are `NA` in each.
#' @param n_taxa Integer vector, how many taxa the table holds.
#'
#' @return `NULL`, invisibly. Errors when `cols` is non-empty.
#'
#' @keywords internal
#' @noRd
mc_stop_partial_na <- function(cols, n_na, n_taxa) {
    if (length(cols) == 0) {
        return(invisible(NULL))
    }

    detail <- paste0(
        cols, " (", n_na, " of ", n_taxa, " taxa NA)"
    )

    stop(
        length(cols), " sample column(s) are only partly measured: ",
        mc_fmt_some(detail),
        ". A sample must be either fully observed or entirely NA, because ",
        "the method imputes whole missing samples rather than scattered ",
        "cells. Correct the data first: set the whole column to NA if that ",
        "sample is genuinely missing, or fill in the individual NA values ",
        "if it is not.",
        call. = FALSE
    )
}

#' Check the table for partly measured samples
#'
#' @description
#' The same check as [mc_stop_partial_na()], applied to a table that has not
#' been through the missing-sample scan, so that `mc_prepare()` refuses the
#' data as well as `mc_detect_missing()`.
#'
#' @param dat The user's data frame.
#' @param col_map Data frame of parsed columns, with a `col` column.
#'
#' @return `NULL`, invisibly. Errors on partly measured columns.
#'
#' @keywords internal
#' @noRd
mc_check_partial_na <- function(dat, col_map) {
    keep <- mc_informative_taxa(dat, col_map)
    n_taxa <- sum(keep)
    if (n_taxa == 0) {
        return(invisible(NULL))
    }

    n_na <- vapply(
        dat[keep, col_map$col, drop = FALSE],
        function(x) sum(is.na(x)), numeric(1)
    )
    bad <- n_na > 0 & n_na < n_taxa

    mc_stop_partial_na(
        col_map$col[bad], n_na[bad], rep(n_taxa, sum(bad))
    )
    invisible(NULL)
}

#' Which taxa carry any data at all
#'
#' @description
#' A taxon that is `NA` everywhere explains its own missing values, so its row
#' must be excluded before judging whether a *column* is only partly measured.
#' Without this, a single empty taxon would make every sample in the table
#' look partly measured.
#'
#' @param dat The user's data frame.
#' @param col_map Data frame of parsed columns, with a `col` column.
#'
#' @return Logical vector, one entry per row of `dat`.
#'
#' @keywords internal
#' @noRd
mc_informative_taxa <- function(dat, col_map) {
    observed <- vapply(
        dat[col_map$col], function(x) !is.na(x), logical(nrow(dat))
    )
    if (is.null(dim(observed))) {
        observed <- matrix(observed, nrow = nrow(dat))
    }
    rowSums(observed) > 0
}

#' Warn about cells that could not be imputed
#'
#' @description
#' Shared by [mc_fit()] and the [mc_run()] reporting path so that a direct
#' call to `mc_fit()` reports unimputable cells too, rather than returning a
#' column of `NA` values in silence.
#'
#' @param pred Long prediction table carrying an `imputed_value` column.
#'
#' @return Integer count of cells left as `NA`, invisibly.
#'
#' @keywords internal
#' @noRd
mc_warn_unimputed <- function(pred) {
    n_failed <- sum(is.na(pred$imputed_value))
    if (n_failed > 0) {
        warning(
            n_failed, " of ", nrow(pred),
            " cell(s) could not be imputed and remain ",
            "NA. This happens when too few subjects have enough ",
            "observations for FPCA to fit that taxon.",
            call. = FALSE
        )
    }
    invisible(n_failed)
}

#' Refuse a mask that leaves nothing to learn from
#'
#' @param mask_pairs Data frame of masked pairs that survived filtering.
#' @param col_map Data frame of parsed columns.
#'
#' @return `NULL`, invisibly. Errors when every sample is masked.
#'
#' @keywords internal
#' @noRd
mc_check_mask_coverage <- function(mask_pairs, col_map) {
    if (nrow(mask_pairs) >= nrow(col_map)) {
        stop(
            "The mask marks every one of the ", nrow(col_map),
            " samples as missing, leaving no observed data to impute from. ",
            "Mask a subset instead.",
            call. = FALSE
        )
    }
    invisible(NULL)
}

#' Refuse arguments the method does not use
#'
#' @description
#' The table methods take `...` so the generic can dispatch, but they do not
#' forward it. Anything extra would be dropped in silence, which is how a
#' misspelled argument, or one renamed in a later version, turns into a run
#' that quietly did something else.
#'
#' `K` is named explicitly because it was the argument for the number of
#' clusters and is now `C`. Silently ignoring it would leave a caller with
#' clustering they did not ask for.
#'
#' @param ... The unused arguments, as received by the method.
#'
#' @return `NULL`, invisibly, when nothing was passed.
#'
#' @keywords internal
#' @noRd
mc_check_dots <- function(...) {
    extra <- names(list(...))
    if (length(extra) == 0) {
        return(invisible(NULL))
    }
    extra[!nzchar(extra)] <- "<unnamed>"

    if ("K" %in% extra) {
        stop(
            "`K` is now `C`, the number of clusters. Passing `K` would ",
            "have been ignored, and the run would have chosen the number ",
            "of clusters itself instead of using the value you gave. ",
            "Rename it to `C`.",
            call. = FALSE
        )
    }
    stop(
        "Unused argument(s): ", paste(extra, collapse = ", "),
        ". Check the spelling against ?mc_run.",
        call. = FALSE
    )
}
