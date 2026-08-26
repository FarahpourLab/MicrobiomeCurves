#' Extract imputed values from FPCA model
#'
#' @description
#' Extracts imputed values and associated information from a fitted
#' FPCA-based longitudinal microbiome model.
#'
#' @details
#' This function returns a data frame containing the true values (if available),
#' imputed values, and squared errors for all masked subject-timepoint
#'  combinations.
#'
#' It is typically used after \code{\link{tti_fit}} to retrieve
#' imputation results for downstream analysis or evaluation.
#'
#' @param fit An object of class \code{"tti_fit"} returned by
#'   \code{\link{tti_fit}}.
#'
#' @return A data frame with columns:
#' \itemize{
#'   \item \code{species}: Taxon identifier
#'   \item \code{rep}: Replicate (subject)
#'   \item \code{time}: Time point
#'   \item \code{true_value}: True value before masking
#'   \item \code{imputed_value}: Imputed value from FPCA
#'   \item \code{FPCA_used}: Logical indicating whether FPCA was successfully
#'   applied
#'   \item \code{se}: Squared error of the imputation
#' }
#'
#' @examples
#' data(taxa_demo)
#'
#' prep <- tti_prepare(
#'     dat = taxa_demo,
#'     taxon_col = "OTU_ID",
#'     mask_list = data.frame(rep = "S01", time = 3)
#' )
#'
#' fit <- suppressWarnings(tti_fit(prep, K = 1))
#'
#' imputed <- tti_impute(fit)
#'
#' head(imputed[, c("species", "rep", "time", "true_value", "imputed_value")])
#'
#' @seealso
#' \code{\link{tti_fit}},
#' \code{\link{tti_metrics}}
#'
#' @export
tti_impute <- function(fit) {

    if (!inherits(fit, "tti_fit")) {
        stop("Input must be an object returned by tti_fit()")
    }

    fit$pred_long
}
