############################################
# One taxon's trajectory, observed and imputed together.
#
# The completed table holds observed and imputed values side by side with
# nothing to tell them apart, and mc_uncertainty() reports only the imputed
# cells. Neither answers the question actually asked of a taxon: what did
# this taxon do over time in each subject, which of those numbers were
# measured, and how well determined were the rest.
#
# This assembles that from what the run already carries. Nothing is refitted.
############################################

#' A taxon's trajectory over time, observed and imputed, with intervals
#'
#' @description
#' Returns one row per subject and time point for the taxa asked for: the
#' value, whether it was observed or imputed, and for an imputed value the
#' 95% interval around it.
#'
#' @details
#' `source` is read from the run's own record of what it imputed, not from
#' whether a column had to be created. The difference matters: a sample whose
#' column was present but entirely `NA` was imputed too, and
#' `run$metadata$imputed` is `FALSE` for it, so a flag taken from there would
#' miss exactly those cells.
#'
#' Why a cell was missing is not reported. `run$design$missing` already
#' says so per cell, and `sample` implies it: a name from your own metadata
#' means the column was there and held `NA`, while a name the run created
#' means no sample existed.
#'
#' The interval is the one the fit produced, read back from the run rather
#' than recomputed, so it agrees exactly with `imputed_uncertainty.tsv`. It
#' is a prediction interval: it includes the residual variance, so it is
#' wider than a confidence interval on the fitted curve, which is the right
#' width for asking how much to trust a single imputed number.
#'
#' The standard error is not reported. A run stores the bounds, not the
#' standard error behind them, and recovering one from the bounds would not
#' reproduce it exactly. [mc_uncertainty()] refits the taxon and gives `se`
#' for the imputed cells.
#'
#' Observed rows carry `NA` in `lower` and `upper`. They were measured, so
#' the run has no interval to offer for them.
#'
#' On a run split by `group_col`, no group column is added: the merged run
#' records subjects and times, not which group an observed sample sat in.
#' Join the result to your own metadata on `subject` when you need it.
#'
#' @param run An object returned by [mc_run()].
#' @param taxon Character vector of taxon names, as they appear in the first
#'   column of `run$completed`. One name or several.
#'
#' @return A data frame ordered by taxon, subject and time, with columns
#'   `taxon`, `subject`, `sample`, `time`, `time_label`, `value`, `source`,
#'   `lower` and `upper`. `sample` is always the column the value
#'   sits in within `run$completed`; for a subject-timepoint that had no
#'   sample of its own that is a name the run created, and `run$metadata`
#'   flags those with `imputed = TRUE`.
#'
#' @examples
#' # A corner of the cohort, so the example is quick.
#' demo <- mc_demo_data(n_taxa = 3, n_subjects = 6)
#' run <- suppressWarnings(mc_run(
#'     demo$abundance, demo$metadata,
#'     sample_col = "sample", subject_col = "subject", time_col = "time",
#'     C = 1, verbose = FALSE
#' ))
#'
#' tr <- mc_imputed_trajectory(run, rownames(demo$abundance)[1])
#' head(tr)
#'
#' # just the imputed cells, widest interval first
#' imp <- tr[tr$source == "imputed", ]
#' imp[order(imp$lower - imp$upper), c("subject", "time", "value",
#'                                     "lower", "upper")]
#'
#' # several taxa at once
#' nrow(mc_imputed_trajectory(run, rownames(demo$abundance)[1:3]))
#'
#' @seealso [mc_uncertainty()] for the imputed cells of one taxon alone,
#'   [mc_run()] for the run that produces these.
#'
#' @export
mc_imputed_trajectory <- function(run, taxon) {
    if (!inherits(run, "mc_run")) {
        stop("run must be an object returned by mc_run().", call. = FALSE)
    }
    if (is.null(run$completed) || ncol(run$completed) < 2L) {
        stop("This run carries no completed table.", call. = FALSE)
    }
    if (is.null(run$metadata) || is.null(run$metadata$subject)) {
        stop(
            "This run carries no per-sample metadata, so subjects and ",
            "times cannot be reported. mc_imputed_trajectory() needs a run ",
            "started from a metadata table.",
            call. = FALSE
        )
    }
    taxon <- mc_trajectory_taxa(run, taxon)

    out <- do.call(rbind, lapply(taxon, function(tx) {
        mc_trajectory_one(run, tx)
    }))
    rownames(out) <- NULL
    out
}


#' Check the taxa asked for against the completed table
#'
#' @param run The `mc_run` object.
#' @param taxon As supplied by the caller.
#'
#' @return `taxon`, as a character vector.
#'
#' @keywords internal
#' @noRd
mc_trajectory_taxa <- function(run, taxon) {
    if (missing(taxon) || is.null(taxon) || length(taxon) == 0) {
        stop(
            "taxon must name at least one taxon. ",
            "Available: ", mc_fmt_some(run$completed[[1]]), ".",
            call. = FALSE
        )
    }
    taxon <- as.character(taxon)
    gone <- setdiff(taxon, run$completed[[1]])
    if (length(gone) > 0) {
        stop(
            "No taxon called ", mc_fmt_some(gone),
            " in the completed table. Available: ",
            mc_fmt_some(run$completed[[1]]), ".",
            call. = FALSE
        )
    }
    dup <- unique(taxon[duplicated(taxon)])
    if (length(dup) > 0) {
        stop(
            "taxon names a taxon more than once: ", mc_fmt_some(dup), ".",
            call. = FALSE
        )
    }
    taxon
}


#' Assemble one taxon's rows
#'
#' @param run The `mc_run` object.
#' @param tx The taxon's name.
#'
#' @return A data frame, one row per subject-timepoint.
#'
#' @keywords internal
#' @noRd
mc_trajectory_one <- function(run, tx) {
    md <- run$metadata
    i <- match(tx, run$completed[[1]])

    # The completed table is the single source for the values themselves, so
    # an observed number here is the caller's own and an imputed one is
    # exactly what was written to disk.
    cols <- match(md$sample, names(run$completed))
    value <- rep(NA_real_, nrow(md))
    have <- !is.na(cols)
    value[have] <- as.numeric(run$completed[i, cols[have]])

    out <- data.frame(
        taxon = tx,
        subject = as.character(md$subject),
        sample = as.character(md$sample),
        time = md$time,
        time_label = if (is.null(md$time_label)) {
            as.character(md$time)
        } else {
            as.character(md$time_label)
        },
        value = value,
        source = "observed",
        lower = NA_real_,
        upper = NA_real_,
        stringsAsFactors = FALSE
    )

    out <- mc_trajectory_mark(out, run, tx)
    out[order(
        match(out$subject, unique(md$subject)), out$time,
        method = "radix"
    ), , drop = FALSE]
}


#' Mark the imputed cells and attach their intervals
#'
#' @description
#' Which cells were imputed is taken from the run's long table, not from
#' whether a column had to be created: a column that existed but held `NA`
#' for every taxon was imputed as well, and is flagged `FALSE` in
#' `run$metadata$imputed`.
#'
#' @param out The taxon's rows so far.
#' @param run The `mc_run` object.
#' @param tx The taxon's name.
#'
#' @return `out`, with `source`, `lower` and `upper` filled in.
#'
#' @keywords internal
#' @noRd
mc_trajectory_mark <- function(out, run, tx) {
    imp <- run$imputed
    if (is.null(imp) || nrow(imp) == 0) {
        return(out)
    }
    imp <- as.data.frame(imp, stringsAsFactors = FALSE)
    imp <- imp[imp$species == tx, , drop = FALSE]
    if (nrow(imp) == 0) {
        return(out)
    }

    key <- function(subject, time) paste(subject, time, sep = "\r")
    j <- match(
        key(out$subject, out$time),
        key(as.character(imp$subject), imp$time)
    )
    hit <- !is.na(j)
    if (!any(hit)) {
        return(out)
    }

    out$source[hit] <- "imputed"
    # Only the bounds. run$imputed carries an se column, but it is a
    # leftover that is never filled, so reading it would promise a number
    # that is always missing.
    for (nm in c("lower", "upper")) {
        if (!is.null(imp[[nm]])) {
            out[[nm]][hit] <- as.numeric(imp[[nm]])[j[hit]]
        }
    }
    out
}
