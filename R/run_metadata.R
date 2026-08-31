# Running from an abundance table plus metadata.
#
# The caller never sees the internal column encoding. Their sample names go
# in, the design is reported, the fit runs on the encoded table, and the
# completed matrix comes back with their sample names restored. Columns
# created for a subject-timepoint that had no sample cannot carry a name of
# theirs, so one is built from the subject and the time and recorded in the
# returned metadata with imputed = TRUE.
#
# Warnings raised along the way are collected rather than left to scroll
# past, so they can be written into the run log beside the design.

#' Impute from an abundance table and its metadata
#'
#' @param dat The abundance table.
#' @param metadata The metadata data.frame.
#' @param sample_col,subject_col,time_col Column names within `metadata`.
#' @param abundance_type Either `"clr"` or `"raw"`.
#' @param pseudocount Zero replacement used when transforming raw values.
#' @param out_dir Directory to write the log and completed table into, or
#'   `NULL` to write nothing.
#' @param K,cluster_method,use_outliers,seed Passed to the fit.
#' @param min_observed Integer. Subjects with fewer observed time points
#'   stop the run.
#' @param verbose Logical. Whether to report progress.
#'
#' @return An object of class `mc_run`, with `completed` carrying the
#'   caller's sample names.
#'
#' @keywords internal
#' @noRd
mc_run_from_metadata <- function(dat, metadata, sample_col, subject_col,
                                  time_col, abundance_type, pseudocount,
                                  out_dir, plots, dpi, plot_format,
                                  K, cluster_method, use_outliers,
                                  seed, min_observed, verbose) {
    mc_check_meta_args(sample_col, subject_col, time_col)
    mc_check_dpi(dpi)

    seen <- new.env(parent = emptyenv())
    seen$warned <- character(0)

    run <- withCallingHandlers(
        {
            design <- mc_from_metadata(
                abundance = dat, metadata = metadata,
                sample_col = sample_col, subject_col = subject_col,
                time_col = time_col, abundance_type = abundance_type,
                pseudocount = pseudocount, verbose = verbose
            )

            fitted <- mc_run_wide(
                design$table,
                taxon_col = design$taxon_col,
                K = K, cluster_method = cluster_method,
                use_outliers = use_outliers, seed = seed,
                min_observed = min_observed, verbose = verbose,
                subject_label = function(codes) {
                    design$subjects[as.integer(sub("^s", "", codes))]
                }
            )

            mc_restore_names(fitted, design)
        },
        warning = function(w) {
            seen$warned <- c(seen$warned, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )

    warned <- seen$warned
    run$warnings <- warned
    mc_replay_warnings(warned)

    if (!is.null(out_dir)) {
        run <- mc_emit_output(run, out_dir, warned, plots, dpi,
                               plot_format, verbose)
    }
    run
}

#' Write everything a finished run has to offer
#'
#' @param run The `mc_run` object.
#' @param out_dir Directory to write into.
#' @param warned Character vector of warnings raised during the run.
#' @param plots Logical. Whether to draw per-taxon uncertainty.
#' @param dpi Resolution for PNG output.
#' @param plot_format One of `"pdf"`, `"png"` or `"both"`.
#' @param verbose Logical. Whether to report what was written.
#'
#' @return `run`, with `files` listing everything written.
#'
#' @keywords internal
#' @noRd
mc_emit_output <- function(run, out_dir, warned, plots, dpi, plot_format,
                            verbose) {
    say <- function(...) if (isTRUE(verbose)) message(...)

    files <- mc_write_output(run, out_dir, warned)
    if (isTRUE(plots)) {
        files <- c(
            files,
            mc_write_uncertainty(run, out_dir, dpi, plot_format, say)
        )
    }

    run$files <- files
    say("Wrote ", length(files), " file(s) to ", out_dir, ": ",
        mc_fmt_some(basename(files), n = 6))
    run
}

#' Re-raise the warnings that were collected for the log
#'
#' @description
#' They were muffled so they could be captured; the caller still needs to
#' see them, so they are raised again once collection is finished.
#'
#' @param warned Character vector of warning messages.
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
#' @noRd
mc_replay_warnings <- function(warned) {
    for (w in warned) {
        warning(w, call. = FALSE)
    }
    invisible(NULL)
}

#' Require all three metadata column arguments together
#'
#' @param sample_col,subject_col,time_col The arguments as supplied.
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
#' @noRd
mc_check_meta_args <- function(sample_col, subject_col, time_col) {
    given <- c(
        sample_col = !missing(sample_col) && !is.null(sample_col),
        subject_col = !missing(subject_col) && !is.null(subject_col),
        time_col = !missing(time_col) && !is.null(time_col)
    )
    if (all(given)) {
        return(invisible(NULL))
    }
    stop(
        "sample_col, subject_col and time_col must all be named, so the ",
        "metadata columns can be identified. Missing: ",
        mc_fmt_some(names(given)[!given]), ".",
        call. = FALSE
    )
}

#' Put the caller's sample names back on a completed table
#'
#' @param run The `mc_run` object produced on the encoded table.
#' @param design The `mc_design` used to encode it.
#'
#' @return The `mc_run` object, renamed, with `design` added.
#'
#' @keywords internal
#' @noRd
mc_restore_names <- function(run, design) {
    completed <- run$completed
    cols <- setdiff(names(completed), design$taxon_col)

    known <- match(cols, design$map$column)
    labels <- ifelse(
        is.na(known),
        mc_name_added(cols, design),
        design$map$sample[known]
    )

    names(completed) <- c("taxon", labels)

    run$completed <- mc_order_samples(completed, cols, design)
    run$imputed <- mc_restore_long(run$imputed, design)
    run$design <- design
    run$metadata <- mc_extend_metadata(design, cols[is.na(known)], labels)
    run$missing <- design$missing
    run$observed <- design$observed
    run
}

#' Put the caller's names back on the long prediction table
#'
#' @description
#' The fit works in encoded subject codes and time positions. Handing those
#' back would make the long table unreadable next to the caller's own
#' metadata, so both are translated and the sample name is added.
#'
#' @param pred The long table from the fit.
#' @param design The `mc_design` the run was built from.
#'
#' @return `pred`, with `subject`, `time` and `sample` in the caller's terms.
#'
#' @keywords internal
#' @noRd
mc_restore_long <- function(pred, design) {
    if (is.null(pred) || nrow(pred) == 0) {
        return(pred)
    }

    pred$subject <- design$subjects[as.integer(sub("^s", "", pred$rep))]
    pred$time <- design$times[pred$time + 1L]
    pred$time_label <- mc_time_label(design$axis, pred$time)

    key <- paste(pred$subject, pred$time, sep = "\r")
    known <- match(key, paste(design$map$subject, design$map$time, sep = "\r"))
    pred$sample <- design$map$sample[known]

    pred$rep <- NULL
    front <- c("species", "subject", "sample", "time", "time_label")
    pred[, c(front, setdiff(names(pred), front)), drop = FALSE]
}

#' Put the completed table into subject and time order
#'
#' @description
#' A column created for an absent sample would otherwise land at the far
#' right, away from the rest of its subject. Ordering by subject and then
#' time puts every sample where a reader expects to find it.
#'
#' @param completed The completed table, already renamed.
#' @param cols Character vector of internal column names, in table order.
#' @param design The `mc_design` the run was built from.
#'
#' @return `completed`, with its sample columns reordered.
#'
#' @keywords internal
#' @noRd
mc_order_samples <- function(completed, cols, design) {
    parsed <- mc_parse_cols(cols)
    subject <- design$subjects[as.integer(sub("^s", "", parsed$rep))]
    time <- design$times[parsed$time + 1L]

    ord <- order(match(subject, design$subjects), time, method = "radix")
    completed[
        , c("taxon", setdiff(names(completed), "taxon")[ord]),
        drop = FALSE
    ]
}
