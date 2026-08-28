# Aggregate reporting for FPCA fits.
#
# mc_safe_fpca() is called once per cluster per masked cell, so on a real
# table it runs tens of thousands of times. Reporting a failure at the call
# site would bury the console, and fdapace's own per-call notes ("The sample
# size is less or equal to 3 curves") repeat just as often. Both are counted
# here instead and reported once, by mc_fit(), when the run finishes.

mc_diag <- new.env(parent = emptyenv())

#' Start a fresh diagnostic tally
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
#' @noRd
mc_diag_reset <- function() {
    mc_diag$calls <- 0L
    mc_diag$failed <- 0L
    mc_diag$small <- 0L
    mc_diag$too_few <- 0L
    mc_diag$kfd_fallback <- 0L
    mc_diag$notes <- list()
    invisible(NULL)
}

#' Record a note raised by the fitting engine
#'
#' @description
#' fdapace raises the same handful of warnings once per fit. Tallying them by
#' text lets [mc_diag_report()] say what was raised and how often, in place
#' of thousands of identical lines.
#'
#' @param text The warning's message.
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
#' @noRd
mc_diag_note <- function(text) {
    if (is.null(mc_diag$calls)) {
        return(invisible(NULL))
    }
    key <- trimws(gsub("[[:space:]]+", " ", text))
    seen <- mc_diag$notes[[key]]
    mc_diag$notes[[key]] <- if (is.null(seen)) 1L else seen + 1L
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
mc_diag_bump <- function(what) {
    if (is.null(mc_diag[[what]])) {
        return(invisible(NULL))
    }
    mc_diag[[what]] <- mc_diag[[what]] + 1L
    invisible(NULL)
}

#' Read the current tally
#'
#' @return A named list of counts, or `NULL` when no tally is active.
#'
#' @keywords internal
#' @noRd
mc_diag_get <- function() {
    if (is.null(mc_diag$calls)) {
        return(NULL)
    }
    list(
        calls = mc_diag$calls,
        failed = mc_diag$failed,
        small = mc_diag$small,
        too_few = mc_diag$too_few,
        kfd_fallback = mc_diag$kfd_fallback,
        notes = mc_diag$notes
    )
}

#' Clear the tally so counting stops
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
#' @noRd
mc_diag_clear <- function() {
    mc_diag$calls <- NULL
    mc_diag$failed <- NULL
    mc_diag$small <- NULL
    mc_diag$too_few <- NULL
    mc_diag$kfd_fallback <- NULL
    mc_diag$notes <- NULL
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
mc_diag_report <- function(n_cells = NA_integer_) {
    d <- mc_diag_get()
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

    parts <- c(parts, mc_diag_note_parts(d$notes))

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

#' Summarise the notes the fitting engine raised
#'
#' @param notes Named list of counts, keyed by message text.
#'
#' @return Character vector, one entry per distinct note, or empty.
#'
#' @keywords internal
#' @noRd
mc_diag_note_parts <- function(notes) {
    if (is.null(notes) || length(notes) == 0) {
        return(character(0))
    }

    counts <- unlist(notes, use.names = TRUE)
    counts <- sort(counts, decreasing = TRUE)

    vapply(
        seq_along(counts),
        function(i) {
            paste0(
                counts[[i]], " fit(s) reported \"",
                mc_truncate(names(counts)[i], 90), "\""
            )
        },
        character(1)
    )
}

#' Shorten a string for use in a message
#'
#' @param x A single string.
#' @param n Maximum number of characters to keep.
#'
#' @return `x`, elided with an ellipsis when longer than `n`.
#'
#' @keywords internal
#' @noRd
mc_truncate <- function(x, n) {
    if (nchar(x) <= n) {
        return(x)
    }
    paste0(substr(x, 1, n - 3), "...")
}
