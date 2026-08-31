# Per-taxon uncertainty of the imputed values.
#
# An imputed number on its own invites more confidence than it deserves. For
# each taxon this computes the analytic interval around every value that was
# imputed, and draws them together, so a taxon whose gaps rest on very little
# information is visible at a glance rather than buried in the table.

#' Interval around every imputed value of one taxon
#'
#' @param fit The `mc_fit` underlying a run.
#' @param species_name Character name of the taxon.
#' @param design The `mc_design` the run was built from, used to report
#'   subjects and times as the caller wrote them.
#'
#' @return Data frame with `subject`, `time`, `time_label`, `imputed`,
#'   `lower`, `upper` and `se`, one row per imputed cell. Rows whose interval
#'   could not be computed carry `NA` bounds.
#'
#' @keywords internal
#' @noRd
mc_taxon_uncertainty <- function(fit, species_name, design) {
    cells <- fit$pred_long[fit$pred_long$species == species_name, ]
    if (nrow(cells) == 0) {
        return(NULL)
    }

    clusters <- fit$clusters[[species_name]]
    W_sp <- mc_subject_time_matrix(fit, species_name)

    bounds <- lapply(seq_len(nrow(cells)), function(i) {
        mc_cell_interval(
            fit, W_sp, clusters, cells$rep[i], cells$time[i]
        )
    })
    bounds <- do.call(rbind, bounds)

    subject <- design$subjects[as.integer(sub("^s", "", cells$rep))]
    time <- design$times[cells$time + 1L]

    data.frame(
        subject = subject,
        time = time,
        time_label = mc_time_label(design$axis, time),
        imputed = cells$imputed_value,
        lower = bounds[, "lower"],
        upper = bounds[, "upper"],
        se = bounds[, "se"],
        stringsAsFactors = FALSE
    )
}

#' Analytic interval for one imputed cell
#'
#' @param fit The `mc_fit`.
#' @param W_sp Subject-by-time matrix for the taxon.
#' @param clusters Named cluster assignment for the taxon, or `NULL`.
#' @param rep_id Character subject code.
#' @param time_id Numeric time position.
#'
#' @return A one-row matrix with `lower`, `upper` and `se`, all `NA` when the
#'   interval cannot be produced.
#'
#' @keywords internal
#' @noRd
mc_cell_interval <- function(fit, W_sp, clusters, rep_id, time_id) {
    none <- cbind(lower = NA_real_, upper = NA_real_, se = NA_real_)
    if (is.null(clusters)) {
        return(none)
    }

    # A subject screened out as an outlier does not appear in the stored
    # clustering, but it was still imputed: the fit keeps the target subject
    # even when flagged, falling back to the full set. The interval follows
    # the same rule, rather than reporting nothing for exactly the subjects
    # whose values are least certain.
    members <- if (rep_id %in% names(clusters)) {
        names(clusters)[clusters == clusters[rep_id]]
    } else {
        union(names(clusters), rep_id)
    }
    fpca <- tryCatch(
        mc_cluster_fpca(W_sp, members, rep_id, time_id, fit$times),
        error = function(e) NULL
    )
    if (is.null(fpca) || is.null(fpca$fp)) {
        return(none)
    }

    ci <- mc_analytic_ci(
        fpca$fp,
        obs_times = fpca$Lt[[rep_id]],
        obs_values = fpca$Ly[[rep_id]],
        pred_time = time_id,
        include_noise = TRUE
    )
    if (is.null(ci)) {
        return(none)
    }
    cbind(lower = ci$lower, upper = ci$upper, se = ci$se)
}

#' Intervals around the imputed values of one taxon
#'
#' @description
#' The numbers behind the uncertainty pages: every value imputed for a taxon
#' with its 95% analytic interval. Useful for filtering a completed table by
#' how well determined each imputed value actually was.
#'
#' @param run An object returned by [mc_run()].
#' @param taxon Character name of the taxon, as it appears in the row names
#'   of the abundance table.
#'
#' @return A data frame with one row per imputed value of that taxon:
#'   `subject`, `time`, `time_label`, `imputed`, `lower`, `upper` and `se`,
#'   with the subject and time given as the caller wrote them.
#'
#' @examples
#' demo <- mc_demo_data()
#' run <- suppressWarnings(mc_run(
#'     demo$counts, demo$metadata,
#'     sample_col = "sample", subject_col = "subject", time_col = "time",
#'     C = 1, verbose = FALSE
#' ))
#'
#' mc_uncertainty(run, rownames(demo$counts)[1])
#'
#' @seealso [mc_run()], whose `out_dir` writes these as one page per value.
#'
#' @export
mc_uncertainty <- function(run, taxon) {
    if (!inherits(run, "mc_run")) {
        stop("run must be an object returned by mc_run().", call. = FALSE)
    }
    if (is.null(run$design)) {
        stop(
            "This run carries no design, so subjects and times cannot be ",
            "reported as you wrote them.",
            call. = FALSE
        )
    }

    known <- unique(run$fit$pred_long$species)
    if (!(taxon %in% known)) {
        stop(
            "No taxon called '", taxon, "' was imputed. Imputed taxa are: ",
            mc_fmt_some(known),
            call. = FALSE
        )
    }

    mc_taxon_uncertainty(run$fit, taxon, run$design)
}
