# Centred log-ratio transformation.
#
# The model works on CLR-transformed abundances: it assumes values live on a
# real line, and raw counts do not. Studies usually hold raw counts or
# relative abundances, so the transformation is offered here rather than left
# as a step people have to remember and get right.
#
# CLR is undefined at zero, which microbiome tables are full of, so zeros are
# replaced first. The default is the standard multiplicative replacement by a
# small fraction of the smallest observed non-zero value in the sample.

#' Transform an abundance matrix to centred log-ratios
#'
#' @param mat Numeric matrix, taxa in rows and samples in columns.
#' @param pseudocount Either `"auto"`, replacing zeros per sample with a
#'   fraction of that sample's smallest non-zero value, or a single positive
#'   number added to every entry.
#'
#' @return A matrix the same shape as `mat`, holding CLR values. Columns that
#'   are entirely `NA` are left as they are, since those are the samples to
#'   be imputed.
#'
#' @keywords internal
#' @noRd
mc_clr <- function(mat, pseudocount = "auto") {
    if (any(mat < 0, na.rm = TRUE)) {
        stop(
            "abundance holds negative values, so it cannot be counts or ",
            "relative abundances. If it is already CLR-transformed, pass ",
            "abundance_type = \"clr\".",
            call. = FALSE
        )
    }

    out <- mat
    for (j in seq_len(ncol(mat))) {
        x <- mat[, j]
        if (all(is.na(x))) next
        out[, j] <- mc_clr_column(x, pseudocount)
    }
    out
}

#' Transform one sample to centred log-ratios
#'
#' @param x Numeric vector of abundances for one sample.
#' @param pseudocount As in [mc_clr()].
#'
#' @return Numeric vector of CLR values.
#'
#' @keywords internal
#' @noRd
mc_clr_column <- function(x, pseudocount) {
    x <- mc_replace_zeros(x, pseudocount)
    lx <- log(x)
    lx - mean(lx, na.rm = TRUE)
}

#' Replace zeros so the logarithm is defined
#'
#' @param x Numeric vector of abundances for one sample.
#' @param pseudocount As in [mc_clr()].
#'
#' @return `x` with zeros replaced.
#'
#' @keywords internal
#' @noRd
mc_replace_zeros <- function(x, pseudocount) {
    if (!identical(pseudocount, "auto")) {
        return(x + pseudocount)
    }

    pos <- x[!is.na(x) & x > 0]
    if (length(pos) == 0) {
        stop(
            "a sample holds no positive value, so it cannot be ",
            "CLR-transformed.",
            call. = FALSE
        )
    }

    # The usual multiplicative replacement: below anything observed, but on
    # the same scale, so the transformed zeros do not dominate the sample.
    x[!is.na(x) & x == 0] <- min(pos) * 0.65
    x
}

#' Check that a pseudocount argument is usable
#'
#' @param pseudocount As supplied by the caller.
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
#' @noRd
mc_check_pseudocount <- function(pseudocount) {
    if (identical(pseudocount, "auto")) {
        return(invisible(NULL))
    }
    ok <- is.numeric(pseudocount) && length(pseudocount) == 1 &&
        !is.na(pseudocount) && pseudocount > 0
    if (!ok) {
        stop(
            "pseudocount must be \"auto\" or a single positive number.",
            call. = FALSE
        )
    }
    invisible(NULL)
}
