# Building the relative abundance and count tables for a finished run.
#
# What is knowable depends on what the caller supplied. Raw counts carry
# their own library sizes and their own zeros, so both can be honoured.
# CLR values carry neither, so the depth is implied from the data and the
# zeros cannot be recovered. That difference is recorded in the run log
# rather than hidden.

#' Completed table on the relative abundance scale
#'
#' @param completed The completed table, taxa in a column named `taxon`.
#' @param pseudocount Per-sample zero threshold, or `NULL` when unknown.
#'
#' @return A data.frame the same shape as `completed`, proportions summing
#'   to one within each sample.
#'
#' @keywords internal
#' @noRd
mc_completed_ra <- function(completed, pseudocount = NULL) {
    out <- completed
    cols <- setdiff(names(completed), "taxon")

    for (cl in cols) {
        p <- if (is.null(pseudocount)) {
            NULL
        } else if (length(pseudocount) > 1 && cl %in% names(pseudocount)) {
            pseudocount[[cl]]
        } else {
            pseudocount[[1]]
        }
        out[[cl]] <- mc_clr_to_ra(completed[[cl]], p)
    }
    out
}

#' Put the caller's own values back over the reconstructed ones
#'
#' @description
#' Applies to samples that were collected. The back-transform is lossy,
#' because the mass taken off the structural zeros is spread evenly rather
#' than back where it came from. For a sample the caller supplied there is
#' no need to pay that cost.
#'
#' @param tbl The reconstructed table, taxa in a column named `taxon`.
#' @param truth Matrix of known values, taxa in rows, or `NULL`.
#'
#' @return `tbl` with the known columns replaced.
#'
#' @keywords internal
#' @noRd
mc_use_known <- function(tbl, truth) {
    if (is.null(truth) || ncol(truth) == 0) {
        return(tbl)
    }
    shared <- intersect(colnames(truth), names(tbl))
    for (cl in shared) {
        tbl[[cl]] <- as.numeric(truth[match(tbl$taxon, rownames(truth)), cl])
    }
    tbl
}

#' Completed table on the count scale
#'
#' @param ra The completed table on the relative abundance scale.
#' @param depths Named numeric vector of library sizes, one per sample, or
#'   `NULL` to imply a depth from each sample.
#'
#' @return A data.frame the same shape as `ra`, holding counts.
#'
#' @keywords internal
#' @noRd
mc_completed_counts <- function(ra, depths = NULL, default = NULL) {
    out <- ra
    cols <- setdiff(names(ra), "taxon")

    for (cl in cols) {
        d <- if (!is.null(depths) && cl %in% names(depths)) {
            depths[[cl]]
        } else if (!is.null(default)) {
            default
        } else {
            mc_implied_depth(ra[[cl]])
        }
        out[[cl]] <- mc_ra_to_counts(ra[[cl]], d)
    }
    out
}

#' Library sizes and zero threshold taken from the caller's raw table
#'
#' @description
#' Only available when the caller supplied counts. A sample that was never
#' collected has no library size of its own, so it takes the median of the
#' samples that were, which keeps its counts on the same scale as the rest
#' of the study.
#'
#' @param abundance The raw abundance matrix as supplied.
#' @param design The `mc_design`, for the sample names.
#'
#' @return List with `depths`, a named vector over every sample in the
#'   completed table, and `pseudocount`, the smallest zero threshold used.
#'
#' @keywords internal
#' @noRd
mc_raw_scale_info <- function(abundance, design) {
    mat <- as.matrix(as.data.frame(abundance))
    depths <- colSums(mat, na.rm = TRUE)
    depths[colSums(!is.na(mat)) == 0] <- NA_real_

    # The CLR step replaced each zero with 0.65 of the smallest positive
    # value in its own sample, then took proportions of the total after
    # that replacement. So on the proportion scale a replaced zero sits at
    # 0.65 of where the rarest real taxon sits. The threshold goes between
    # the two, which separates them exactly rather than nearly.
    pseudo <- vapply(seq_len(ncol(mat)), function(j) {
        x <- mat[, j]
        pos <- x[!is.na(x) & x > 0]
        if (length(pos) == 0 || is.na(depths[j]) || depths[j] <= 0) {
            return(NA_real_)
        }
        n_zero <- sum(!is.na(x) & x == 0)
        total <- depths[j] + n_zero * min(pos) * 0.65
        (min(pos) * 0.825) / total
    }, numeric(1))
    names(pseudo) <- colnames(mat)

    # A sample that was collected is kept as it was given. Its counts and
    # its proportions are known, so reconstructing them through the CLR and
    # back would replace the caller's own numbers with worse ones.
    observed <- mat[, !is.na(depths) & depths > 0, drop = FALSE]

    typical <- stats::median(depths, na.rm = TRUE)
    list(
        observed = observed,
        depths = depths,
        default_depth = typical,
        pseudocount = pseudo,
        default_pseudocount = stats::median(pseudo, na.rm = TRUE)
    )
}
