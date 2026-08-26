#' Assemble the data frames the panels draw
#'
#' @param fit Object of class `tti_fit`.
#' @param species_name Character name of the taxon.
#' @param rep_id Character identifier of the target subject.
#' @param time_id Numeric time point being imputed.
#' @param W_sp Subject-by-time matrix.
#' @param fpca List returned by [tti_cluster_fpca()].
#' @param members Character vector of subjects in the target's cluster.
#'
#' @return List of data frames: `fit_cluster`, `target`, `obs`, `mask`, `imp`
#'   and `true`.
#'
#' @keywords internal
#' @noRd
tti_plot_frames <- function(fit, species_name, rep_id, time_id,
                            W_sp, fpca, members) {
    times <- fit$times
    col_map <- fit$col_map
    taxon_col <- fit$taxon_col

    df_fit <- as.data.frame(fpca$fit_mat)
    colnames(df_fit) <- fpca$grid
    df_fit <- dplyr::mutate(df_fit, replicate = names(fpca$Ly))
    df_fit <- tidyr::pivot_longer(
        df_fit, -replicate,
        names_to = "time", values_to = "value"
    )
    df_fit <- dplyr::mutate(df_fit, time = as.numeric(time))

    df_fit_cluster <- dplyr::filter(df_fit, replicate %in% members)
    df_target <- dplyr::filter(df_fit_cluster, replicate == rep_id)

    pts <- tti_plot_points(
        fit, species_name, rep_id, time_id, W_sp, fpca
    )

    list(
        fit_cluster = df_fit_cluster, target = df_target, obs = pts$obs,
        mask = pts$mask, imp = pts$imp, true = pts$true
    )
}

#' Observed, masked, imputed and true points for the target subject
#'
#' @param fit Object of class `tti_fit`.
#' @param species_name Character name of the taxon.
#' @param rep_id Character identifier of the target subject.
#' @param time_id Numeric time point being imputed.
#' @param W_sp Subject-by-time matrix.
#' @param fpca List returned by [tti_cluster_fpca()].
#'
#' @return List of data frames: `obs`, `mask`, `imp` and `true`.
#'
#' @keywords internal
#' @noRd
tti_plot_points <- function(fit, species_name, rep_id, time_id,
                            W_sp, fpca) {
    times <- fit$times
    col_map <- fit$col_map
    taxon_col <- fit$taxon_col

    df_obs <- dplyr::filter(
        tibble::tibble(
            time = times,
            value = W_sp[rep_id, ],
            replicate = rep_id
        ),
        time != time_id
    )

    df_mask <- tibble::tibble(
        time = time_id,
        value = W_sp[rep_id, which(times == time_id)]
    )

    subj_i <- match(rep_id, names(fpca$Ly))
    t_pos <- which.min(abs(fpca$grid - time_id))

    true_val <- dplyr::pull(
        dplyr::filter(fit$dat_orig, .data[[taxon_col]] == species_name),
        col_map$col[col_map$rep == rep_id & col_map$time == time_id]
    )

    df_imp <- tibble::tibble(
        time = time_id, value = fpca$fit_mat[subj_i, t_pos]
    )
    df_true <- tibble::tibble(time = time_id, value = true_val)

    df_obs$type <- "Observed"
    df_imp$type <- "Imputed"
    df_true$type <- "True"

    list(obs = df_obs, mask = df_mask, imp = df_imp, true = df_true)
}

#' Draw the fitted panel with its uncertainty band
#'
#' @param frames List of data frames from [tti_plot_frames()].
#' @param bands List with `analytic` and `bootstrap` interval data frames,
#'   either of which may be `NULL`.
#' @param rep_id Character identifier of the target subject.
#' @param time_id Numeric time point being imputed.
#' @param species_name Character name of the taxon.
#' @param cl_id Cluster identifier of the target subject.
#' @param ylim Numeric length-2 vector giving the y range.
#'
#' @return A `ggplot` object.
#'
#' @keywords internal
#' @noRd
tti_plot_fitted_panel <- function(frames, bands, rep_id, time_id,
                                  species_name, cl_id, ylim) {
    tti_add_ci_ribbons(ggplot(), bands) +
        geom_line(
            data = dplyr::filter(frames$fit_cluster, replicate != rep_id),
            aes(time, value, group = replicate),
            color = "#bdbdbd", linewidth = 1
        ) +
        geom_line(
            data = frames$target, aes(time, value),
            color = "#1f78b4", linewidth = 1.6, linetype = "dashed"
        ) +
        geom_point(
            data = frames$obs, aes(time, value),
            color = "#1f78b4", size = 2.5
        ) +
        geom_point(
            data = frames$mask, aes(time, value),
            shape = 21, fill = "white", color = "#1f78b4", size = 3
        ) +
        geom_point(
            data = frames$imp, aes(time, value),
            color = "#e31a1c", shape = 17, size = 3.5
        ) +
        geom_point(
            data = frames$true, aes(time, value),
            color = "#33a02c", size = 3.5
        ) +
        geom_vline(
            xintercept = time_id, linetype = "dashed", color = "#636363"
        ) +
        coord_cartesian(ylim = ylim) +
        labs(
            title = paste0("Cluster ", cl_id, " FPCA fit with uncertainty"),
            subtitle = paste0(species_name, " | Replicate: ", rep_id),
            x = "Time",
            y = "Centered Log-Ratio"
        ) +
        scale_fill_manual(
            name = "Uncertainty",
            values = c(
                "Bootstrap CI" = "gray70",
                "Analytic CI" = "#80b1d3"
            )
        ) +
        theme_minimal(base_size = 13)
}

#' Draw the raw trajectory panel
#'
#' @param W_sp Subject-by-time matrix.
#' @param outliers Named logical vector marking outlying subjects.
#' @param df_mask One-row data frame holding the masked point.
#' @param rep_id Character identifier of the target subject.
#' @param time_id Numeric time point being imputed.
#' @param ylim Numeric length-2 vector giving the y range.
#'
#' @return A `ggplot` object.
#'
#' @keywords internal
#' @noRd
tti_plot_raw_panel <- function(W_sp, outliers, df_mask, rep_id,
                               time_id, ylim) {
    df_raw <- tti_raw_long(W_sp, outliers)

    df_others <- dplyr::filter(df_raw, replicate != rep_id)
    df_raw_target <- dplyr::filter(df_raw, replicate == rep_id)

    ggplot() +
        geom_line(
            data = df_others,
            aes(time, value,
                group = interaction(replicate),
                color = outlier_flag
            )
        ) +
        geom_point(
            data = df_others,
            aes(time, value, color = outlier_flag)
        ) +
        geom_line(
            data = df_raw_target, aes(time, value),
            color = "#1f78b4", linewidth = 1.6, linetype = "dashed"
        ) +
        geom_point(
            data = df_raw_target, aes(time, value),
            color = "#1f78b4", size = 2.5
        ) +
        geom_point(
            data = df_mask, aes(time, value),
            shape = 21, fill = "white", color = "#1f78b4", size = 3
        ) +
        geom_vline(
            xintercept = time_id, linetype = "dashed", color = "#636363"
        ) +
        coord_cartesian(ylim = ylim) +
        labs(
            title = "Raw trajectories",
            x = "Time",
            y = "Centered Log-Ratio"
        ) +
        scale_color_manual(
            values = c("FALSE" = "#bdbdbd", "TRUE" = "red"),
            name = "Trajectory type",
            labels = c("FALSE" = "Non-outlier", "TRUE" = "Outlier")
        ) +
        theme_minimal(base_size = 13)
}

#' Y ranges for the two panels
#'
#' @param W_sp Subject-by-time matrix.
#' @param members Character vector of subjects in the target's cluster.
#' @param outliers Named logical vector marking outlying subjects.
#' @param fpca List returned by [tti_cluster_fpca()].
#' @param ci_obj List with `analytic` and `bootstrap` data frames.
#'
#' @return List with `left` and `right`, each a length-2 range.
#'
#' @keywords internal
#' @noRd
tti_plot_ylims <- function(W_sp, members, outliers, fpca, ci_obj) {
    left <- range(
        W_sp[members, ],
        fpca$fit_mat,
        if (!is.null(ci_obj$bootstrap)) {
            c(ci_obj$bootstrap$lower, ci_obj$bootstrap$upper)
        },
        if (!is.null(ci_obj$analytic)) {
            c(ci_obj$analytic$lower, ci_obj$analytic$upper)
        },
        na.rm = TRUE
    )

    reps_right <- unique(c(members, names(outliers)[outliers]))
    right <- range(W_sp[reps_right, , drop = FALSE], na.rm = TRUE)

    list(left = left, right = right)
}

#' Long form of the raw trajectories, with outlier flags
#'
#' @param W_sp Subject-by-time matrix.
#' @param outliers Named logical vector marking outlying subjects.
#'
#' @return Data frame with `replicate`, `time`, `value` and
#'   `outlier_flag`.
#'
#' @keywords internal
#' @noRd
tti_raw_long <- function(W_sp, outliers) {
    df_raw <- as.data.frame(W_sp)
    df_raw$replicate <- rownames(df_raw)
    df_raw <- tidyr::pivot_longer(
        df_raw, -replicate,
        names_to = "time", values_to = "value"
    )
    df_raw <- dplyr::mutate(df_raw, time = as.numeric(time))
    df_raw$outlier_flag <- outliers[df_raw$replicate]
    df_raw
}

#' Add the uncertainty ribbons to the fitted panel
#'
#' @param p A `ggplot` object.
#' @param bands List with `analytic` and `bootstrap` data frames.
#'
#' @return The `ggplot` object with any available ribbons added.
#'
#' @keywords internal
#' @noRd
tti_add_ci_ribbons <- function(p, bands) {
    if (!is.null(bands$bootstrap)) {
        boot <- bands$bootstrap
        boot$type <- "Bootstrap CI"
        p <- p + geom_ribbon(
            data = boot,
            aes(x = time, ymin = lower, ymax = upper, fill = type),
            alpha = 0.25
        )
    }

    if (!is.null(bands$analytic)) {
        ana <- bands$analytic
        ana$type <- "Analytic CI"
        p <- p + geom_ribbon(
            data = ana,
            aes(x = time, ymin = lower, ymax = upper, fill = type),
            alpha = 0.35
        )
    }

    p
}
