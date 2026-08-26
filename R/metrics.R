#' @importFrom dplyr group_by summarise
NULL
#' Compute imputation performance metrics
#'
#' @description
#' Computes mean squared error (MSE) and root mean squared error (RMSE)
#' for FPCA-based imputation of longitudinal microbiome data.
#'
#' @details
#' Metrics are calculated by comparing imputed values to the true values
#' at masked subject-timepoint locations.
#'
#' The function provides:
#' \itemize{
#'   \item Overall MSE across all imputed points
#'   \item Overall RMSE
#'   \item MSE for each (replicate, time) combination
#' }
#'
#' These metrics are useful for evaluating imputation accuracy in
#' simulation studies or cross-validation settings where true values
#' are known.
#'
#' @param fit An object of class \code{"tti_fit"} returned by
#'   \code{\link{tti_fit}}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{MSE_overall}: Mean squared error across all imputed values
#'   \item \code{RMSE_overall}: Root mean squared error
#'   \item \code{MSE_by_point}: Data frame with MSE for each (replicate, time)
#' }
#'
#' @examples
#' data(taxa_demo)
#'
#' # Metrics are only meaningful when the true values are known, i.e. after
#' # masking values on purpose. They are NA for a tti_run() fit, where the
#' # missing values were never observed.
#' prep <- tti_prepare(
#'     dat = taxa_demo,
#'     taxon_col = "OTU_ID",
#'     mask_list = data.frame(rep = "S01", time = 3)
#' )
#'
#' fit <- suppressWarnings(tti_fit(prep, K = 1))
#'
#' metrics <- tti_metrics(fit)
#'
#' metrics$MSE_overall
#' metrics$RMSE_overall
#' head(metrics$MSE_by_point)
#'
#' @seealso
#' \code{\link{tti_fit}}
#'
#' @export
tti_metrics <- function(fit) {

    if (!inherits(fit, "tti_fit")) {
        stop("Input must be an object returned by tti_fit()")
    }

    pred_long <- fit$pred_long

    # A fit from tti_run() imputes samples that were absent from the table,
    # so there is no held-out truth to score against and every error is NA.
    # Without this the caller silently receives NaN.
    if (all(is.na(pred_long$true_value))) {
        warning(
            "No true values are available, so every metric is NaN. This fit ",
            "imputes samples that were missing from the data, which have no ",
            "known value to compare against. Accuracy can only be measured ",
            "on a fit whose masked cells were observed, as produced by ",
            "tti_prepare() with mask_list or mask_matrix.",
            call. = FALSE
        )
    }

    MSE_by_point <- dplyr::group_by(pred_long, rep, time)

    MSE_by_point <- dplyr::summarise(
        MSE_by_point,
        MSE = mean(se, na.rm = TRUE),
        .groups = "drop"
    )

    overall_MSE <- mean(pred_long$se, na.rm = TRUE)
    overall_RMSE <- sqrt(overall_MSE)

    list(
        MSE_overall = overall_MSE,
        RMSE_overall = overall_RMSE,
        MSE_by_point = MSE_by_point
    )
}
