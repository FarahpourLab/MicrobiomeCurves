# Per-taxon uncertainty of the imputed values.
#
# An imputed number on its own invites more confidence than it deserves. For
# each taxon this computes the analytic interval around every value that was
# imputed, and draws them together, so a taxon whose gaps rest on very little
# information is visible at a glance rather than buried in the table.

#' Interval around every imputed value of one taxon
#'
#' @param fit The `tti_fit` underlying a run.
#' @param species_name Character name of the taxon.
#' @param design The `tti_design` the run was built from, used to report
#'   subjects and times as the caller wrote them.
#'
#' @return Data frame with `subject`, `time`, `time_label`, `imputed`,
#'   `lower`, `upper` and `se`, one row per imputed cell. Rows whose interval
#'   could not be computed carry `NA` bounds.
#'
#' @keywords internal
#' @noRd
tti_taxon_uncertainty <- function(fit, species_name, design) {
    cells <- fit$pred_long[fit$pred_long$species == species_name, ]
    if (nrow(cells) == 0) {
        return(NULL)
    }

    clusters <- fit$clusters[[species_name]]
    W_sp <- tti_subject_time_matrix(fit, species_name)

    bounds <- lapply(seq_len(nrow(cells)), function(i) {
        tti_cell_interval(
            fit, W_sp, clusters, cells$rep[i], cells$time[i]
        )
    })
    bounds <- do.call(rbind, bounds)

    subject <- design$subjects[as.integer(sub("^s", "", cells$rep))]
    time <- design$times[cells$time + 1L]

    data.frame(
        subject = subject,
        time = time,
        time_label = tti_time_label(design$axis, time),
        imputed = cells$imputed_value,
        lower = bounds[, "lower"],
        upper = bounds[, "upper"],
        se = bounds[, "se"],
        stringsAsFactors = FALSE
    )
}

#' Analytic interval for one imputed cell
#'
#' @param fit The `tti_fit`.
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
tti_cell_interval <- function(fit, W_sp, clusters, rep_id, time_id) {
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
        tti_cluster_fpca(W_sp, members, rep_id, time_id, fit$times),
        error = function(e) NULL
    )
    if (is.null(fpca) || is.null(fpca$fp)) {
        return(none)
    }

    ci <- tti_analytic_ci(
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
