# The trajectory view of one imputed value.
#
# A point with an error bar says how wide the interval is but not why. This
# draws the thing the interval came from: every subject's trajectory, the
# fitted curve through them, and the 95% band around the prediction for the
# subject whose value was missing.
#
# When outlier screening is on, two fits are drawn -- one over all subjects
# and one with the flagged trajectories removed -- because the gap between
# them is what screening actually did to the answer. With screening off
# there is one fit and no flagged trajectories, and the panel says so
# instead of showing an empty distinction.

TTI_COL <- c(
    all = "#D55E00", kept = "#0072B2", retained = "#334155",
    flagged = "#B3261E"
)

#' Draw the trajectories and interval behind one imputed value
#'
#' @param panel List from [tti_panel_data()].
#' @param use_outliers Logical. Whether the fit screened for outliers.
#'
#' @return A `ggplot` object.
#'
#' @keywords internal
#' @noRd
tti_uncertainty_panel <- function(panel, use_outliers) {
    p <- ggplot2::ggplot() +
        ggplot2::geom_ribbon(
            data = panel$bands,
            ggplot2::aes(
                .data$time, ymin = .data$lower, ymax = .data$upper,
                fill = .data$set
            ),
            alpha = 0.13
        ) +
        ggplot2::geom_line(
            data = panel$traj,
            ggplot2::aes(
                .data$time, .data$value,
                group = .data$subject, colour = .data$grp,
                linewidth = .data$grp
            )
        ) +
        ggplot2::geom_line(
            data = panel$bands,
            ggplot2::aes(.data$time, .data$fit, colour = .data$set),
            linewidth = 0.8
        ) +
        ggplot2::geom_vline(
            xintercept = panel$target_time,
            linewidth = 0.3, linetype = "22", colour = "grey40"
        )

    tti_panel_points(p, panel) + tti_panel_scales(panel, use_outliers)
}

#' Add the observed and imputed points to a panel
#'
#' @param p The plot so far.
#' @param panel List from [tti_panel_data()].
#'
#' @return The plot with its points added.
#'
#' @keywords internal
#' @noRd
tti_panel_points <- function(p, panel) {
    p <- p +
        ggplot2::geom_point(
            data = panel$observed,
            ggplot2::aes(.data$time, .data$value, shape = "Observed"),
            size = 1.3
        ) +
        ggplot2::geom_point(
            data = panel$imputed,
            ggplot2::aes(.data$time, .data$value, shape = "Imputed"),
            size = 2.3, colour = "black", fill = "white", stroke = 0.6
        )

    if (nrow(panel$truth) > 0) {
        # Only a benchmarking fit knows what was masked; a real run does not.
        p <- p + ggplot2::geom_point(
            data = panel$truth,
            ggplot2::aes(.data$time, .data$value, shape = "True value"),
            size = 2.3, colour = "black", fill = "grey70", stroke = 0.6
        )
    }
    p
}

#' Colour, fill and shape scales for a panel
#'
#' @param panel List from [tti_panel_data()].
#' @param use_outliers Logical. Whether the fit screened for outliers.
#'
#' @return A list of ggplot scales and theming.
#'
#' @keywords internal
#' @noRd
tti_panel_scales <- function(panel, use_outliers) {
    cols <- c(
        "Fit, all replicates" = unname(TTI_COL["all"]),
        "Fit, outliers excluded" = unname(TTI_COL["kept"]),
        "Retained" = unname(TTI_COL["retained"]),
        "Flagged outlier" = unname(TTI_COL["flagged"]),
        "Subjects" = unname(TTI_COL["retained"]),
        "Fit" = unname(TTI_COL["kept"])
    )

    shapes <- c("Observed" = 16, "Imputed" = 21, "True value" = 21)

    list(
        ggplot2::scale_colour_manual(
            values = cols, breaks = panel$colour_breaks, name = NULL
        ),
        ggplot2::scale_fill_manual(
            values = cols, breaks = panel$fill_breaks, name = "95% interval"
        ),
        ggplot2::scale_shape_manual(values = shapes, name = NULL),
        ggplot2::scale_linewidth_manual(
            values = c(
                "Retained" = 0.45, "Flagged outlier" = 0.8,
                "Subjects" = 0.45
            ),
            guide = "none"
        ),
        ggplot2::scale_x_continuous(
            breaks = panel$axis_breaks, labels = panel$axis_labels
        ),
        ggplot2::labs(
            title = panel$title, subtitle = panel$subtitle,
            x = "Time", y = "CLR abundance"
        ),
        ggplot2::theme_classic(base_size = 9),
        ggplot2::theme(
            legend.position = "bottom", legend.box = "vertical",
            plot.title = ggplot2::element_text(size = 9, hjust = 0),
            plot.subtitle = ggplot2::element_text(size = 8, hjust = 0)
        )
    )
}
