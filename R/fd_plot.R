#' @importFrom ggplot2 ggplot geom_line geom_point geom_ribbon geom_vline
#' @importFrom ggplot2 coord_cartesian labs scale_fill_manual scale_color_manual
#'  theme_minimal aes
#' @importFrom ggplot2 coord_cartesian labs scale_fill_manual theme_minimal aes
#' @importFrom dplyr filter mutate pull
#' @importFrom patchwork wrap_plots
#' @importFrom tidyr pivot_longer
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
NULL
#' Plot FPCA-based imputation with uncertainty
#'
#' @description
#' Visualizes FPCA-based imputation for longitudinal microbiome data,
#' including fitted trajectories, observed values, imputed values, and
#' uncertainty bands.
#'
#' @details
#' This function generates a two-panel plot:
#'
#' \itemize{
#'   \item \strong{Left panel}: FPCA-fitted trajectories with uncertainty bands
#'   \item \strong{Right panel}: Raw observed trajectories
#' }
#'
#' The left panel includes:
#' \itemize{
#'   \item Cluster member trajectories (gray lines)
#'   \item Target subject trajectory (blue dashed line)
#'   \item Observed data points (blue)
#'   \item Masked point (open circle)
#'   \item Imputed value (red triangle)
#'   \item True value (green point, if available)
#'   \item Confidence bands:
#'     \itemize{
#'       \item Analytic (model-based)
#'       \item Bootstrap (resampling-based)
#'     }
#' }
#'
#' The FPCA model is fitted using only subjects within the same cluster
#' as the target subject.
#'
#' @param fit An object of class \code{"mc_fit"} returned by
#'   \code{\link{mc_fit}}.
#'
#' @param species_name Character. Name of the taxon (feature) to plot.
#'
#' @param rep_id Character. Target subject (replicate).
#'
#' @param time_id Numeric. Time point corresponding to the masked observation.
#'
#' @param ci_method Character. Type of uncertainty to display:
#'   \code{"analytic"}, \code{"bootstrap"}, or \code{"both"}.
#'
#' @param B Integer. Number of bootstrap samples used when applicable.
#'
#' @return A \code{patchwork} object combining FPCA fit and raw data plots.
#'
#' @examples
#' data(taxa_demo)
#'
#' prep <- mc_prepare(
#'     dat = taxa_demo,
#'     taxon_col = "OTU_ID",
#'     mask_list = data.frame(rep = "S01", time = 3)
#' )
#' fit <- suppressWarnings(mc_fit(prep, C = 1))
#'
#' suppressWarnings(
#'     mc_plot(
#'         fit = fit,
#'         species_name = fit$pred_long$species[1],
#'         rep_id = "S01",
#'         time_id = 3,
#'         ci_method = "analytic"
#'     )
#' )
#'
#' # Bootstrap bands refit the model B times. B is small here to keep the
#' # example quick; use a few hundred for real work.
#' suppressWarnings(
#'     mc_plot(
#'         fit = fit,
#'         species_name = fit$pred_long$species[1],
#'         rep_id = "S01",
#'         time_id = 3,
#'         ci_method = "both",
#'         B = 25
#'     )
#' )
#'
#' @seealso
#' \code{\link{mc_fit}},
#' \code{\link{mc_ci}}
#' @export
mc_plot <- function(
    fit,
    species_name,
    rep_id,
    time_id,
    ci_method = c("analytic", "bootstrap", "both"),
    B = 200
) {
    ci_method <- match.arg(ci_method)

    clusters <- fit$clusters[[species_name]]
    outliers <- fit$outliers[[species_name]]

    if (is.null(clusters)) {
        stop("No clustering found for this species")
    }

    cl_id <- clusters[rep_id]
    members <- names(clusters)[clusters == cl_id]

    W_sp <- mc_subject_time_matrix(fit, species_name)
    fpca <- mc_cluster_fpca(W_sp, members, rep_id, time_id, fit$times)

    ci_obj <- mc_ci(
        fit = fit,
        species_name = species_name,
        rep_id = rep_id,
        time_id = time_id,
        method = ci_method,
        B = B
    )

    frames <- mc_plot_frames(
        fit, species_name, rep_id, time_id, W_sp, fpca, members
    )

    ylim <- mc_plot_ylims(W_sp, members, outliers, fpca, ci_obj)

    p1 <- mc_plot_fitted_panel(
        frames, ci_obj, rep_id, time_id, species_name, cl_id, ylim$left
    )
    p2 <- mc_plot_raw_panel(
        W_sp, outliers, frames$mask, rep_id, time_id, ylim$right
    )

    p1 + p2
}
