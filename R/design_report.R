# What the user is told about their design, on screen.
#
# The report is deliberately concrete: counts first, then the specific
# subject-timepoints that are missing, then anything that will weaken the
# fit. Everything printed here is also returned in the tti_design object, so
# nothing is only available by reading the console.

#' Report a design as it is built
#'
#' @param design An object of class `tti_design`.
#' @param say Function used to emit each line.
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
#' @noRd
tti_report_design <- function(design, say) {
    n_samples <- nrow(design$map)
    n_cells <- length(design$subjects) * length(design$times)

    say(
        "Design: ", nrow(design$table), " taxa, ", n_samples, " samples, ",
        length(design$subjects), " subjects, ",
        length(design$times), " time points."
    )
    say("  time points, in order: ", tti_axis_text(design$axis))

    tti_report_gaps(design, n_cells, say)
    tti_report_thin(design, say)
    invisible(NULL)
}

#' Report the missing part of a design
#'
#' @param design An object of class `tti_design`.
#' @param n_cells Integer size of the full subject-by-time grid.
#' @param say Function used to emit each line.
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
#' @noRd
tti_report_gaps <- function(design, n_cells, say) {
    miss <- design$missing
    if (nrow(miss) == 0) {
        say("  complete: every subject has a sample at every time point.")
        return(invisible(NULL))
    }

    say(
        "  missing: ", nrow(miss), " of ", n_cells,
        " subject-timepoints (",
        format(round(100 * nrow(miss) / n_cells, 1), nsmall = 1), "%)."
    )

    for (r in c("absent_sample", "no_data")) {
        sel <- miss[miss$reason == r, , drop = FALSE]
        if (nrow(sel) == 0) next
        say(
            "    ", nrow(sel), " ", tti_reason_text(r), ": ",
            tti_fmt_some(paste0(sel$subject, " at ", sel$time_label), n = 6)
        )
    }
    invisible(NULL)
}

#' Describe a missingness reason in words
#'
#' @param reason One of `"absent_sample"` or `"no_data"`.
#'
#' @return A single string.
#'
#' @keywords internal
#' @noRd
tti_reason_text <- function(reason) {
    switch(reason,
        absent_sample = "with no sample at all",
        no_data = "whose sample column is entirely NA",
        reason
    )
}

#' Warn about subjects that carry little information
#'
#' @param design An object of class `tti_design`.
#' @param say Function used to emit each line.
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
#' @noRd
tti_report_thin <- function(design, say) {
    obs <- design$observed
    thin <- obs$subject[obs$n_observed < 2]
    if (length(thin) == 0) {
        return(invisible(NULL))
    }

    say(
        "  note: ", length(thin),
        " subject(s) have fewer than two observed time points (",
        tti_fmt_some(thin, n = 6),
        "). Their values are drawn largely from the population curve."
    )
    invisible(NULL)
}

#' Print a design
#'
#' @param x An object of class `tti_design`.
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#'
#' @examples
#' counts <- matrix(
#'     c(1.2, 0.4, 1.5, 0.6, 1.1, 0.5),
#'     nrow = 2,
#'     dimnames = list(
#'         c("Bacteroides", "Prevotella"),
#'         c("RUN_0031", "RUN_0044", "RUN_0052")
#'     )
#' )
#' meta <- data.frame(
#'     library = c("RUN_0031", "RUN_0044", "RUN_0052"),
#'     animal = c("M01", "M01", "M02"),
#'     day = c(0, 7, 0)
#' )
#' design <- tti_from_metadata(
#'     counts, meta,
#'     sample_col = "library", subject_col = "animal", time_col = "day",
#'     verbose = FALSE
#' )
#' design
#'
#' @export
print.tti_design <- function(x, ...) {
    cat("TaxaTimeImpute design\n")
    cat("  taxa        :", nrow(x$table), "\n")
    cat("  samples     :", nrow(x$map), "\n")
    cat("  subjects    :", length(x$subjects), "\n")
    cat("  time points :", tti_axis_text(x$axis), "\n")
    cat(
        "  missing     :", nrow(x$missing), "of",
        length(x$subjects) * length(x$times), "subject-timepoints\n"
    )
    invisible(x)
}
