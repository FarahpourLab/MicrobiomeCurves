# Assembling everything one uncertainty panel needs.
#
# The band is the prediction interval for the subject whose value was
# missing, evaluated across the fitted grid rather than at the single
# imputed point, so the panel shows how the certainty changes either side of
# the gap.

#' Gather the pieces behind one imputed value
#'
#' @param fit The `tti_fit` underlying a run.
#' @param species_name Character name of the taxon.
#' @param rep_id Character subject code, internal.
#' @param time_id Numeric time position.
#' @param design The `tti_design`, for the caller's own names.
#' @param use_outliers Logical. Whether the fit screened for outliers.
#'
#' @return A list of data frames and labels for [tti_uncertainty_panel()], or
#'   `NULL` when no band could be produced.
#'
#' @keywords internal
#' @noRd
tti_panel_data <- function(fit, species_name, rep_id, time_id, design,
                           use_outliers) {
    W_sp <- tti_subject_time_matrix(fit, species_name)
    outliers <- fit$outliers[[species_name]]
    flagged <- if (isTRUE(use_outliers) && !is.null(outliers)) {
        names(outliers)[outliers]
    } else {
        character(0)
    }

    bands <- tti_panel_bands(fit, W_sp, rep_id, time_id, flagged)
    if (is.null(bands)) {
        return(NULL)
    }

    subj <- function(x) design$subjects[as.integer(sub("^s", "", x))]
    target <- subj(rep_id)
    tlab <- tti_time_label(design$axis, design$times[time_id + 1L])

    list(
        bands = bands,
        traj = tti_panel_traj(W_sp, design, flagged, use_outliers),
        observed = tti_panel_observed(W_sp, rep_id, design),
        imputed = tti_panel_imputed(fit, species_name, rep_id, time_id,
                                    design),
        truth = tti_panel_truth(fit, species_name, rep_id, time_id, design),
        target_time = time_id,
        cell_label = paste0(target, " at ", tlab),
        axis_breaks = seq_along(design$times) - 1L,
        axis_labels = tti_time_label(design$axis, design$times),
        title = paste0(species_name, ": ", target, " at time ", tlab),
        subtitle = tti_panel_subtitle(flagged, subj, use_outliers),
        colour_breaks = if (length(flagged) > 0) {
            c("Fit, all replicates", "Fit, outliers excluded",
              "Retained", "Flagged outlier")
        } else if (isTRUE(use_outliers)) {
            c("Fit, all replicates", "Retained")
        } else {
            c("Fit", "Subjects")
        },
        fill_breaks = unique(bands$set)
    )
}

#' Say in words what screening did to this taxon
#'
#' @param flagged Character vector of flagged subject codes.
#' @param subj Function mapping codes to the caller's subject names.
#' @param use_outliers Logical. Whether screening was on.
#'
#' @return A single string.
#'
#' @keywords internal
#' @noRd
tti_panel_subtitle <- function(flagged, subj, use_outliers) {
    if (!isTRUE(use_outliers)) {
        return("outlier screening off; one fit over all subjects")
    }
    if (length(flagged) == 0) {
        return("outlier screening on; no trajectory was flagged")
    }
    paste0(
        "outlier screening on; flagged: ",
        tti_fmt_some(subj(flagged), n = 5)
    )
}
