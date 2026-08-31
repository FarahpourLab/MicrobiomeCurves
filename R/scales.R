# Returning the completed table on the scales people work with.
#
# The model works in centred log ratios, so that is what it produces. Most
# downstream analysis wants relative abundance, and some tools want counts.
# Both are derived here rather than left to the user, because the
# back-transform has a step that is easy to get wrong: the inverse of a CLR
# is strictly positive, so every structural zero comes back as a small
# positive number unless it is put back deliberately.
#
# The route follows the evaluation code used for the published benchmark:
# invert the CLR, zero anything below the pseudo-count, then spread the mass
# that was removed across the taxa that remain.

#' Back-transform centred log ratios to relative abundance
#'
#' @param clr Numeric vector for one sample.
#' @param pseudocount Value below which a back-transformed proportion is
#'   treated as a structural zero, or `NULL` to keep every value positive.
#'
#' @return Numeric vector of proportions summing to one.
#'
#' @keywords internal
#' @noRd
mc_clr_to_ra <- function(clr, pseudocount = NULL) {
    if (all(is.na(clr))) {
        return(clr)
    }

    # Inverse CLR is the softmax: exp(x) shifted for numerical safety, then
    # scaled to sum to one.
    e <- exp(clr - max(clr, na.rm = TRUE))
    ra <- e / sum(e, na.rm = TRUE)

    if (is.null(pseudocount)) {
        return(ra)
    }

    ra[!is.na(ra) & ra < pseudocount] <- 0
    nz <- !is.na(ra) & ra > 0
    if (any(nz)) {
        # The zeroed mass is redistributed evenly over the taxa that remain,
        # so the sample still sums to one.
        ra[nz] <- ra[nz] + (1 - sum(ra, na.rm = TRUE)) / sum(nz)
    }
    ra
}

#' Library size implied by a relative abundance vector
#'
#' @description
#' The rarest taxon that is present must have been seen at least once, so
#' the reciprocal of the smallest non-zero proportion is a lower bound on
#' the number of reads. Capped, because a very small proportion implies an
#' implausible depth.
#'
#' @param ra Numeric vector of proportions.
#' @param cap Integer upper bound on the implied depth.
#'
#' @return A single integer, or `NA_integer_` when nothing is positive.
#'
#' @keywords internal
#' @noRd
mc_implied_depth <- function(ra, cap = 100000L) {
    pos <- ra[!is.na(ra) & ra > 0]
    if (length(pos) == 0) {
        return(NA_integer_)
    }
    as.integer(min(cap, round(1 / min(pos))))
}

#' Turn relative abundance into counts at a given depth
#'
#' @description
#' Rounded rather than sampled. The evaluation code draws a multinomial,
#' which suits simulating diversity, but a results file that changes
#' between runs of the same data would be surprising.
#'
#' @param ra Numeric vector of proportions.
#' @param depth Library size to scale to.
#'
#' @return Numeric vector of counts.
#'
#' @keywords internal
#' @noRd
mc_ra_to_counts <- function(ra, depth) {
    if (is.na(depth) || depth <= 0) {
        return(rep(NA_real_, length(ra)))
    }
    round(ra * depth)
}
