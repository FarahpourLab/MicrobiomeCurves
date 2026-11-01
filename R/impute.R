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
#' It is typically used after \code{\link{mc_fit}} to retrieve
#' imputation results for downstream analysis or evaluation.
#'
#' @param fit An object of class \code{"mc_fit"} returned by
#'   \code{\link{mc_fit}}.
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
#' data(taxa_demo_meta)
#'
#' # A corner of the cohort, so the example is quick: three taxa and six
#' # subjects, of the 833 and 64 the bundled table holds.
#' subs <- head(unique(taxa_demo_meta$subject), 6)
#' md <- taxa_demo_meta[taxa_demo_meta$subject %in% subs, ]
#' small <- taxa_demo[1:3, c("OTU_ID", md$sample)]
#'
#' # A subject and time the subset has, to mask and then recover.
#' cell <- data.frame(rep = md$subject[1], time = md$time[1])
#'
#' prep <- mc_prepare(
#'     dat = small,
#'     taxon_col = "OTU_ID",
#'     mask_list = cell
#' )
#'
#' fit <- suppressWarnings(mc_fit(prep, C = 1))
#'
#' imputed <- mc_impute(fit)
#'
#' head(imputed[, c("species", "rep", "time", "true_value", "imputed_value")])
#'
#' @seealso
#' \code{\link{mc_fit}},
#' \code{\link{mc_metrics}}
#'
#' @export
mc_impute <- function(fit) {

    if (!inherits(fit, "mc_fit")) {
        stop("Input must be an object returned by mc_fit()")
    }

    fit$pred_long
}
