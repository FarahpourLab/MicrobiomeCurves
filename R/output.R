# Writing a run to disk.
#
# What is written is the record of the run: the design, the time points in
# order, every subject-timepoint that was missing and why, any warning
# raised, and the completed table itself. Someone opening the folder months
# later should be able to see what was imputed and on what basis, without
# rerunning anything.

#' Write a run's log and completed table to a directory
#'
#' @param run The `tti_run` object.
#' @param out_dir Directory to write into. Created if it does not exist.
#' @param warnings Character vector of warnings raised during the run.
#'
#' @return Character vector of the two paths written, invisibly.
#'
#' @keywords internal
#' @noRd
tti_write_output <- function(run, out_dir, warnings = character(0)) {
    if (!dir.exists(out_dir)) {
        ok <- dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        if (!ok) {
            stop(
                "Could not create out_dir '", out_dir, "'.",
                call. = FALSE
            )
        }
    }

    tsv <- file.path(out_dir, "imputed_abundance.tsv")
    log <- file.path(out_dir, "imputation_log.txt")

    utils::write.table(
        run$completed, tsv,
        sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
    )
    writeLines(tti_log_lines(run, warnings), log)

    invisible(c(tsv, log))
}

#' Compose the lines of the run log
#'
#' @param run The `tti_run` object.
#' @param warnings Character vector of warnings raised during the run.
#'
#' @return Character vector of lines.
#'
#' @keywords internal
#' @noRd
tti_log_lines <- function(run, warnings) {
    d <- run$design

    c(
        "TaxaTimeImpute run log",
        paste0("package version: ", tti_version_string()),
        "",
        "DESIGN",
        paste0("  taxa        : ", nrow(run$completed)),
        paste0("  samples     : ", nrow(d$map)),
        paste0("  subjects    : ", length(d$subjects)),
        paste0("  time points : ", length(d$times)),
        paste0("  time order  : ", tti_axis_text(d$axis)),
        paste0("  subjects    : ", tti_fmt_some(d$subjects, n = 20)),
        "",
        tti_log_missing(run, d),
        "",
        tti_log_result(run),
        "",
        tti_log_warnings(warnings)
    )
}

#' The missing-sample section of the log
#'
#' @param run The `tti_run` object.
#' @param d Its design.
#'
#' @return Character vector of lines.
#'
#' @keywords internal
#' @noRd
tti_log_missing <- function(run, d) {
    n_cells <- length(d$subjects) * length(d$times)
    miss <- d$missing

    head <- paste0(
        "MISSING (", nrow(miss), " of ", n_cells, " subject-timepoints)"
    )
    if (nrow(miss) == 0) {
        return(c(head, "  none: every subject has a sample at every time."))
    }

    body <- vapply(
        seq_len(nrow(miss)),
        function(i) {
            paste0(
                "  ", miss$subject[i], "  at time ", miss$time_label[i],
                "  [", miss$reason[i], "]",
                if (!is.na(miss$sample[i])) {
                    paste0("  sample ", miss$sample[i])
                } else {
                    "  no sample was collected"
                }
            )
        },
        character(1)
    )
    c(head, body)
}

#' The outcome section of the log
#'
#' @param run The `tti_run` object.
#'
#' @return Character vector of lines.
#'
#' @keywords internal
#' @noRd
tti_log_result <- function(run) {
    n_cells <- nrow(run$imputed)
    created <- sum(run$metadata$imputed)

    c(
        "RESULT",
        paste0("  cells imputed  : ", n_cells - run$n_failed, " of ", n_cells),
        paste0("  columns added  : ", created),
        paste0(
            "  added columns  : ",
            if (created > 0) {
                tti_fmt_some(run$metadata$sample[run$metadata$imputed], n = 20)
            } else {
                "none"
            }
        )
    )
}

#' The warnings section of the log
#'
#' @param warnings Character vector of warnings raised during the run.
#'
#' @return Character vector of lines.
#'
#' @keywords internal
#' @noRd
tti_log_warnings <- function(warnings) {
    if (length(warnings) == 0) {
        return(c("WARNINGS", "  none"))
    }
    c(
        paste0("WARNINGS (", length(warnings), ")"),
        paste0("  - ", gsub("[[:space:]]+", " ", warnings))
    )
}

#' The package version, for the log header
#'
#' @return A single string.
#'
#' @keywords internal
#' @noRd
tti_version_string <- function() {
    v <- tryCatch(
        as.character(utils::packageVersion("TaxaTimeImpute")),
        error = function(e) "unknown"
    )
    v
}
