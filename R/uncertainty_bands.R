# The pieces of a panel: bands, trajectories and points.

#' Prediction bands across the fitted grid
#'
#' @description
#' One band when screening is off or nothing was flagged, two when it was:
#' the fit over all subjects and the fit with the flagged ones removed. The
#' gap between them is what screening did to this value.
#'
#' @param fit The `mc_fit`.
#' @param W_sp Subject-by-time matrix for the taxon.
#' @param rep_id Character subject code of the target.
#' @param time_id Numeric time position being imputed.
#' @param flagged Character vector of flagged subject codes.
#'
#' @return Data frame with `time`, `fit`, `lower`, `upper` and `set`, or
#'   `NULL` when no band could be produced.
#'
#' @keywords internal
#' @noRd
mc_panel_bands <- function(fit, W_sp, rep_id, time_id, flagged) {
    everyone <- rownames(W_sp)
    kept <- setdiff(everyone, flagged)
    if (!(rep_id %in% kept)) {
        kept <- union(kept, rep_id)
    }

    sets <- if (length(flagged) > 0 && length(kept) >= 2) {
        list(
            "Fit, all replicates" = everyone,
            "Fit, outliers excluded" = kept
        )
    } else {
        stats::setNames(list(everyone), if (length(flagged) > 0) {
            "Fit, all replicates"
        } else {
            "Fit"
        })
    }

    bands <- lapply(names(sets), function(nm) {
        b <- mc_one_band(fit, W_sp, sets[[nm]], rep_id, time_id)
        if (is.null(b)) NULL else cbind(b, set = nm)
    })
    bands <- Filter(Negate(is.null), bands)
    if (length(bands) == 0) {
        return(NULL)
    }
    do.call(rbind, bands)
}

#' One prediction band over the fitted grid
#'
#' @param fit The `mc_fit`.
#' @param W_sp Subject-by-time matrix for the taxon.
#' @param members Character vector of subjects to fit on.
#' @param rep_id Character subject code of the target.
#' @param time_id Numeric time position being imputed.
#'
#' @return Data frame with `time`, `fit`, `lower` and `upper`, or `NULL`.
#'
#' @keywords internal
#' @noRd
mc_one_band <- function(fit, W_sp, members, rep_id, time_id) {
    fpca <- tryCatch(
        mc_cluster_fpca(W_sp, members, rep_id, time_id, fit$times),
        error = function(e) NULL
    )
    if (is.null(fpca) || is.null(fpca$fp)) {
        return(NULL)
    }

    grid <- fpca$fp$workGrid
    ci <- lapply(grid, function(tt) {
        mc_analytic_ci(
            fpca$fp,
            obs_times = fpca$Lt[[rep_id]], obs_values = fpca$Ly[[rep_id]],
            pred_time = tt, include_noise = TRUE
        )
    })
    ok <- !vapply(ci, is.null, logical(1))
    if (!any(ok)) {
        return(NULL)
    }

    data.frame(
        time = grid[ok],
        fit = vapply(ci[ok], `[[`, numeric(1), "mean"),
        lower = vapply(ci[ok], `[[`, numeric(1), "lower"),
        upper = vapply(ci[ok], `[[`, numeric(1), "upper"),
        stringsAsFactors = FALSE
    )
}

#' Every subject's observed trajectory
#'
#' @param W_sp Subject-by-time matrix for the taxon.
#' @param design The `mc_design`.
#' @param flagged Character vector of flagged subject codes.
#' @param use_outliers Logical. Whether screening was on.
#'
#' @return Data frame with `subject`, `time`, `value` and `grp`.
#'
#' @keywords internal
#' @noRd
mc_panel_traj <- function(W_sp, design, flagged, use_outliers) {
    codes <- rownames(W_sp)
    times <- as.numeric(colnames(W_sp))

    out <- lapply(codes, function(s) {
        v <- as.numeric(W_sp[s, ])
        keep <- !is.na(v)
        if (!any(keep)) {
            return(NULL)
        }
        data.frame(
            subject = design$subjects[as.integer(sub("^s", "", s))],
            time = times[keep], value = v[keep],
            grp = if (!isTRUE(use_outliers)) {
                "Subjects"
            } else if (s %in% flagged) {
                "Flagged outlier"
            } else {
                "Retained"
            },
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, Filter(Negate(is.null), out))
}

#' The target subject's own observations
#'
#' @param W_sp Subject-by-time matrix for the taxon.
#' @param rep_id Character subject code of the target.
#' @param design The `mc_design`.
#'
#' @return Data frame with `time` and `value`.
#'
#' @keywords internal
#' @noRd
mc_panel_observed <- function(W_sp, rep_id, design) {
    v <- as.numeric(W_sp[rep_id, ])
    keep <- !is.na(v)
    data.frame(
        time = as.numeric(colnames(W_sp))[keep], value = v[keep],
        stringsAsFactors = FALSE
    )
}

#' The imputed value itself
#'
#' @param fit The `mc_fit`.
#' @param species_name,rep_id,time_id The cell.
#' @param design The `mc_design`.
#'
#' @return One-row data frame with `time` and `value`.
#'
#' @keywords internal
#' @noRd
mc_panel_imputed <- function(fit, species_name, rep_id, time_id, design) {
    pl <- fit$pred_long
    row <- pl$species == species_name & pl$rep == rep_id &
        pl$time == time_id
    data.frame(
        time = time_id,
        value = pl$imputed_value[row][1],
        stringsAsFactors = FALSE
    )
}

#' The masked value, when the fit knows it
#'
#' @param fit The `mc_fit`.
#' @param species_name,rep_id,time_id The cell.
#' @param design The `mc_design`.
#'
#' @return A data frame with no rows when the true value is unknown, which is
#'   the case for any run over genuinely missing data.
#'
#' @keywords internal
#' @noRd
mc_panel_truth <- function(fit, species_name, rep_id, time_id, design) {
    pl <- fit$pred_long
    row <- pl$species == species_name & pl$rep == rep_id &
        pl$time == time_id
    truth <- pl$true_value[row][1]

    if (is.na(truth)) {
        return(data.frame(time = numeric(0), value = numeric(0)))
    }
    data.frame(
        time = time_id, value = truth,
        stringsAsFactors = FALSE
    )
}
