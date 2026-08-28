#' @importFrom dplyr filter distinct left_join select mutate bind_rows
#' @importFrom tidyr pivot_longer
#' @importFrom tibble tibble rownames_to_column
#' @importFrom stringr str_match
NULL
#' Prepare longitudinal microbiome data for FPCA imputation
#'
#' @description
#' Parses and validates longitudinal microbiome data with replicate-time
#'  structure,
#' and defines masked subject-timepoint combinations where entire samples are
#'  missing.
#'
#' This function is the first step of the MicrobiomeCurves workflow. It extracts
#'  replicate
#' and time information from column names, identifies valid observations, and
#'  constructs
#' a structured object used for downstream FPCA-based imputation.
#'
#' @details
#' The function expects column names (except the taxon column) to follow a
#'  pattern such as:
#' \code{replicate.time} (e.g., "A1.1", "A1.2", "A2.1").
#'
#' Users can define missing samples (entire subject-timepoints) using:
#' \itemize{
#'   \item \code{mask_list}: a data frame or list of (replicate, time) pairs
#'   \item \code{mask_matrix}: a matrix indicating observed (1) or missing (0)
#'   samples
#' }
#'
#' These masked entries represent situations where all taxa/features are missing
#'  for a given
#' subject at a specific time point.
#'
#' @param dat A data.frame containing microbiome measurements in wide format.
#'   Rows correspond to taxa/features and columns correspond to replicate-time
#'   observations.
#'
#' @param taxon_col A character string specifying the column containing taxon
#'  identifiers
#'   (e.g., OTU IDs, species names). Default is \code{"OTU_ID"}.
#'
#' @param mask_list Optional. A data.frame or list specifying missing samples as
#'  pairs
#'   of (replicate, time).
#'
#'   As a data frame, the first two columns are read as replicate and time,
#'   whatever they are called:
#'   \code{data.frame(rep = c("A1"), time = c(3))}.
#'
#'   As a list, each element is one pair of length two, so the replicate sits
#'   inside the element: \code{list(c("A1", 3), c("A2", 5))}. Note that
#'   \code{list(A1 = 3)} is a named number rather than a pair and is
#'   rejected.
#'
#' @param mask_matrix Optional. A matrix/data.frame where rows are replicates
#'  and columns
#'   are time points, with values 1 (observed) or 0 (missing).
#'
#' @param parse_fun Optional. A custom function to parse column names into
#'  replicate and time.
#'   It must return a data.frame with columns \code{col}, \code{rep}, and
#'   \code{time}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{dat}: Original input data
#'   \item \code{taxon_col}: Name of taxon column
#'   \item \code{col_map}: Mapping of columns to replicate and time
#'   \item \code{reps}: Unique replicate identifiers
#'   \item \code{times}: Sorted unique time points
#'   \item \code{mask_pairs}: Valid masked replicate-time combinations
#'   \item \code{parse_fun}: Parsing function used
#' }
#'
#' @examples
#' data(taxa_demo)
#'
#' # Benchmark style: hide values that are known, to score them later.
#' mask_list <- data.frame(
#'     rep  = c("S01", "S02"),
#'     time = c(3, 5)
#' )
#'
#' prep <- mc_prepare(
#'     dat = taxa_demo,
#'     taxon_col = "OTU_ID",
#'     mask_list = mask_list
#' )
#'
#' prep$mask_pairs
#' head(prep$col_map)
#'
#' # Use mc_run() to impute samples that are missing in the data. It
#' # detects them and does not need a mask.
#'
#' @export
mc_prepare <- function(
    dat,
    taxon_col = "OTU_ID",
    mask_list = NULL,
    mask_matrix = NULL,
    parse_fun = NULL
) {
    if (is.null(parse_fun)) {
        parse_fun <- mc_parse_cols
    }

    col_map <- mc_prepare_col_map(dat, taxon_col, parse_fun)
    reps <- unique(col_map$rep)
    times <- sort(unique(col_map$time))

    mc_check_inputs(dat, taxon_col, col_map, reps)

    mask_pairs <- dplyr::distinct(
        mc_collect_mask_pairs(mask_list, mask_matrix)
    )
    requested <- mask_pairs
    mask_pairs <- dplyr::filter(mask_pairs, rep %in% reps, time %in% times)
    if (nrow(mask_pairs) == 0) {
        mc_stop_no_mask_pairs(requested, reps, times)
    }
    mc_check_mask_coverage(mask_pairs, col_map)

    mask_pairs <- mask_pairs %>%
        dplyr::left_join(col_map, by = c("rep", "time")) %>%
        dplyr::filter(!is.na(col))

    list(
        dat = dat,
        taxon_col = taxon_col,
        col_map = col_map,
        reps = reps,
        times = times,
        mask_pairs = mask_pairs,
        parse_fun = parse_fun
    )
}
