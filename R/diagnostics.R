# Aggregate reporting for FPCA fits.
#
# tti_safe_fpca() is called once per cluster per masked cell, so on a real
# table it runs tens of thousands of times. Reporting a failure at the call
# site would bury the console, and fdapace's own per-call notes ("The sample
# size is less or equal to 3 curves") repeat just as often. Both are counted
# here instead and reported once, by tti_fit(), when the run finishes.

tti_diag <- new.env(parent = emptyenv())

#' Start a fresh diagnostic tally
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
#' @noRd
tti_diag_reset <- function() {
    tti_diag$calls <- 0L
    tti_diag$failed <- 0L
    tti_diag$small <- 0L
    tti_diag$too_few <- 0L
    tti_diag$kfd_fallback <- 0L
    invisible(NULL)
}

#' Record the outcome of one FPCA attempt
#'
#' @param what One of `"calls"`, `"failed"`, `"small"` or `"too_few"`.
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
#' @noRd
tti_diag_bump <- function(what) {
    if (is.null(tti_diag[[what]])) {
        return(invisible(NULL))
    }
    tti_diag[[what]] <- tti_diag[[what]] + 1L
    invisible(NULL)
}

#' Read the current tally
#'
#' @return A named list of counts, or `NULL` when no tally is active.
#'
#' @keywords internal
#' @noRd
tti_diag_get <- function() {
    if (is.null(tti_diag$calls)) {
        return(NULL)
    }
    list(
        calls = tti_diag$calls,
        failed = tti_diag$failed,
        small = tti_diag$small,
        too_few = tti_diag$too_few,
        kfd_fallback = tti_diag$kfd_fallback
    )
}

#' Clear the tally so counting stops
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
#' @noRd
tti_diag_clear <- function() {
    tti_diag$calls <- NULL
    tti_diag$failed <- NULL
    tti_diag$small <- NULL
    tti_diag$too_few <- NULL
    tti_diag$kfd_fallback <- NULL
    invisible(NULL)
}

#' Report the tally as a single warning
#'
#' @description
#' Says how often the functional fit could not be produced, and how often it
#' rested on very few curves. Both mean the affected cells fall back to a
#' wider subject set or are left unimputed, which the user cannot otherwise
#' see.
#'
#' @param n_cells Integer number of masked cells in the run, used to put the
#'   failure count in context.
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
#' @noRd
tti_diag_report <- function(n_cells = NA_integer_) {
    d <- tti_diag_get()
    if (is.null(d) || d$calls == 0L) {
        return(invisible(NULL))
    }

    parts <- character(0)
    if (d$too_few > 0L) {
        parts <- c(parts, paste0(
            d$too_few, " fit(s) were skipped because fewer than two ",
            "subjects had at least two observations"
        ))
    }
    if (d$failed > 0L) {
        parts <- c(parts, paste0(
            d$failed, " fit(s) failed to converge"
        ))
    }
    if (d$small > 0L) {
        parts <- c(parts, paste0(
            d$small, " fit(s) used three or fewer curves, so the ",
            "functional basis is poorly determined"
        ))
    }
    if (isTRUE(d$kfd_fallback > 0L)) {
        parts <- c(parts, paste0(
            d$kfd_fallback, " functional k-means call(s) failed and fell ",
            "back to clustering the FPC scores, which kmeans.fd() does on ",
            "small sets of curves"
        ))
    }

    if (length(parts) > 0) {
        warning(
            "FPCA diagnostics across ", d$calls, " fit(s): ",
            paste(parts, collapse = "; "),
            ". Affected cells fall back to the full subject set or are ",
            "left as NA; check FPCA_used in the result.",
            call. = FALSE
        )
    }
    invisible(NULL)
}
