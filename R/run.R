############################################
# User-facing pipeline
#
# This file is additive. It does not modify tti_prepare(), tti_fit(),
# tti_impute() or any part of the FPCA / CI computation. It only
#
#   (a) works out which subject-timepoints are missing in a user's table,
#   (b) hands them to the existing functions in exactly the form the
#       benchmark drivers already use, and
#   (c) writes the returned predictions back into the user's layout.
#
# Numerical results are therefore identical to calling
# tti_prepare() -> tti_fit() -> tti_impute() by hand.
############################################


#' Default column-name builder
#'
#' Inverse of \code{\link{tti_parse_cols}}, used to name columns that have to
#' be created for subject-timepoints absent from the input table.
#'
#' @param rep Character vector of subject identifiers.
#' @param time Numeric vector of time points.
#' @return Character vector of column names.
#' @keywords internal
#' @noRd
tti_default_make_col <- function(rep, time) {
    paste0(rep, ".", time)
}


#' Find missing subject-timepoints in a longitudinal microbiome table
#'
#' @description
#' Inspects a wide abundance table and reports which whole samples
#' (subject-timepoint combinations) are missing, without fitting anything.
#'
#' @details
#' Two kinds of missingness are recognised, both of which mean "this entire
#' sample was not observed":
#'
#' \itemize{
#'   \item \code{"absent_column"}  -  the subject-timepoint has no column in the
#'     table at all.
#'   \item \code{"all_na"}  -  the column exists but every taxon in it is
#'     \code{NA}.
#' }
#'
#' The intended sampling grid is taken to be every subject crossed with every
#' time point seen in the column names, unless \code{times} is given
#' explicitly. Inspect the result before calling \code{\link{tti_run}} if you
#' are unsure whether that assumption fits your design.
#'
#' A sample column must be either fully observed or entirely \code{NA}. The
#' method models whole missing samples, so a column holding some taxa and
#' \code{NA} for others is neither, and is rejected with an error naming the
#' offending columns. Correct such a column before calling: make it entirely
#' \code{NA} if the sample is genuinely missing, or fill in the individual
#' \code{NA} values if it is not.
#'
#' @param dat A data.frame in wide format: one column of taxon identifiers,
#'   the rest named \code{"<subject>.<time>"}.
#' @param taxon_col Character. Name of the taxon identifier column.
#' @param times Optional numeric vector giving the intended time grid. Defaults
#'   to the sorted unique times found in the column names.
#' @param parse_fun Optional custom parser, as in \code{\link{tti_prepare}}.
#'   Must return columns \code{col}, \code{rep}, \code{time}.
#' @param make_col Optional function of \code{(rep, time)} returning a column
#'   name, used only when a subject-timepoint has to be created. Defaults to
#'   \code{paste0(rep, ".", time)}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{missing}: data.frame of \code{rep}, \code{time}, \code{reason}
#'   \item \code{observed}: data.frame of \code{rep} and \code{n_observed}
#'   \item \code{col_map}: the parsed column map
#'   \item \code{reps}, \code{times}: the subject and time grids
#'   \item \code{partial_na}: data.frame of columns that are partly \code{NA}.
#'   Always empty in a returned value, since such a column stops the call;
#'   the field is retained so the structure does not depend on the data.
#' }
#'
#' @examples
#' data(taxa_demo)
#'
#' # The bundled example table already contains missing samples.
#' info <- tti_detect_missing(taxa_demo, taxon_col = "OTU_ID")
#'
#' # what is missing, and why
#' info$missing
#'
#' # how many observations each subject has left
#' head(info$observed)
#'
#' # the same table read from a CSV instead
#' path <- system.file("extdata", "taxa_demo.csv", package = "TaxaTimeImpute")
#' dat <- read.csv(path, check.names = FALSE)
#' nrow(tti_detect_missing(dat, taxon_col = "OTU_ID")$missing)
#'
#' @seealso \code{\link{tti_run}}
#' @export
tti_detect_missing <- function(
    dat,
    taxon_col = "OTU_ID",
    times = NULL,
    parse_fun = NULL,
    make_col = NULL
) {
    if (is.null(parse_fun)) parse_fun <- tti_parse_cols
    if (is.null(make_col)) make_col <- tti_default_make_col

    parsed <- tti_column_map(dat, taxon_col, parse_fun)
    col_map <- parsed$col_map

    reps <- unique(col_map$rep)
    times <- if (is.null(times)) {
        sort(unique(col_map$time))
    } else {
        sort(unique(times))
    }

    found <- tti_find_missing(dat, col_map, reps, times, make_col)
    missing <- found$missing

    # These are returned in the result either way, but a direct caller of
    # tti_detect_missing() would otherwise never be told they exist.
    tti_warn_column_issues(parsed$unparsed, found$partial_na)

    observed <- tti_observed_counts(missing, reps, times)

    list(
        missing = missing,
        observed = observed,
        col_map = col_map,
        reps = reps,
        times = times,
        partial_na = found$partial_na,
        unparsed = parsed$unparsed,
        taxon_col = taxon_col
    )
}


#' Impute missing timepoints in a longitudinal microbiome table
#'
#' @description
#' One-call pipeline. Takes a wide abundance table containing missing whole
#' samples, works out what is missing, fits the FPCA model, and returns the
#' table with the missing samples filled in.
#'
#' @details
#' This is a convenience wrapper. It performs no modelling of its own: it
#' calls \code{\link{tti_prepare}}, \code{\link{tti_fit}} and
#' \code{\link{tti_impute}} exactly as the benchmark scripts do, so the numbers
#' it returns are identical to running those three functions by hand.
#'
#' Missingness is detected by \code{\link{tti_detect_missing}} rather than
#' supplied as a mask, which is the difference between this function and the
#' benchmark workflow. Because the true values are unknown, the
#' \code{true_value} and \code{se} columns of the returned long table are
#' \code{NA}, and \code{\link{tti_metrics}} is not meaningful on the fit.
#'
#' \strong{Subjects with few observations.} A subject needs at least
#' \code{min_observed} observed time points for its own trajectory to inform
#' its imputation; below that the prediction falls back towards the population
#' mean curve. Such subjects are still imputed and still contribute to the
#' model fit, but a warning names them so the results can be interpreted
#' accordingly. This mirrors the rule used throughout the published benchmark.
#'
#' @param dat An abundance table. With \code{metadata}, this is the natural
#'   form: taxa in rows, one column per sample, named however the study names
#'   them. Without \code{metadata}, the sample columns must instead be named
#'   \code{"<subject>.<time>"}. Values are expected to be CLR-transformed
#'   abundances.
#' @param taxon_col Character. Name of the taxon identifier column. With
#'   \code{metadata} it may be \code{NULL}, in which case taxa are taken from
#'   the row names.
#' @param metadata Optional data.frame with one row per sample, saying which
#'   subject each sample came from and when. Supplying it lets \code{dat}
#'   keep the study's own sample names, which are also used for the returned
#'   \code{completed} table. The three columns below must then be named.
#' @param sample_col Character. Name of the \code{metadata} column holding
#'   sample identifiers, matching the columns of \code{dat}.
#' @param subject_col Character. Name of the \code{metadata} column
#'   identifying whom or what the repeated measurements were taken on: a
#'   mouse, a participant, a plot, a bioreactor.
#' @param time_col Character. Name of the \code{metadata} column holding the
#'   time point. Values must be numeric or coercible to numeric.
#' @param K Integer or \code{NULL}. Number of trajectory clusters, passed
#'   straight to \code{\link{tti_fit}}. \code{NULL} selects K per taxon by mean
#'   silhouette width.
#' @param cluster_method Character, passed to \code{\link{tti_fit}}.
#' @param use_outliers Logical, passed to \code{\link{tti_fit}}. \code{TRUE}
#'   enables functional-depth outlier screening.
#' @param seed Integer random seed, passed to \code{\link{tti_fit}}.
#' @param times Optional numeric vector giving the intended time grid, passed
#'   to \code{\link{tti_detect_missing}}.
#' @param parse_fun,make_col Optional column-name parser and builder, passed to
#'   \code{\link{tti_detect_missing}} and \code{\link{tti_prepare}}.
#' @param min_observed Integer. Observation count below which a subject is
#'   flagged in a warning. Default \code{2}.
#' @param verbose Logical. Print progress. Default \code{TRUE}.
#'
#' @return An object of class \code{"tti_run"}: a list with
#' \itemize{
#'   \item \code{completed}: the input table with missing samples filled in.
#'     Original columns keep their order and their values; columns created for
#'     absent subject-timepoints are appended.
#'   \item \code{imputed}: long table of every imputed cell (taxon x subject x
#'     time) with \code{imputed_value} and \code{FPCA_used}.
#'   \item \code{missing}: what was detected as missing, and why.
#'   \item \code{observed}: observation count per subject.
#'   \item \code{n_failed}: number of cells FPCA could not impute (left
#'     \code{NA} in \code{completed}).
#'   \item \code{fit}: the underlying \code{"tti_fit"} object.
#' }
#'
#' @examples
#' data(taxa_demo)
#'
#' # A complete run on the bundled example data.
#' # fdapace emits "time gap" notices on sparse designs; see ?taxa_demo.
#' run <- suppressWarnings(
#'     tti_run(taxa_demo, taxon_col = "OTU_ID", K = 1, verbose = FALSE)
#' )
#'
#' run
#'
#' # the filled table, in the same layout as the input
#' run$completed[1:4, 1:5]
#'
#' # every imputed cell, in long form
#' head(run$imputed[, c("species", "rep", "time", "imputed_value")])
#'
#' # what was treated as missing
#' run$missing
#'
#' # A real analysis starts from a file of your own. The bundled CSV stands
#' # in for one here, and the result is written to a temporary directory.
#' path <- system.file("extdata", "taxa_demo.csv", package = "TaxaTimeImpute")
#' dat <- read.csv(path, check.names = FALSE)
#'
#' run2 <- suppressWarnings(
#'     tti_run(dat, taxon_col = "OTU_ID", K = 1, verbose = FALSE)
#' )
#'
#' out <- file.path(tempdir(), "taxa_demo_imputed.csv")
#' write.csv(run2$completed, out, row.names = FALSE)
#' file.exists(out)
#'
#' @seealso
#' \code{\link{tti_detect_missing}},
#' \code{\link{tti_prepare}},
#' \code{\link{tti_fit}},
#' \code{\link{tti_impute}}
#'
#' @export
setGeneric("tti_run", function(dat, ...) standardGeneric("tti_run"))

#' @rdname tti_run
#' @export
setMethod("tti_run", "data.frame", function(
    dat,
    taxon_col = "OTU_ID",
    metadata = NULL,
    sample_col = NULL,
    subject_col = NULL,
    time_col = NULL,
    K = NULL,
    cluster_method = c("fpca", "kmeans_fd"),
    use_outliers = TRUE,
    seed = 123,
    times = NULL,
    parse_fun = NULL,
    make_col = NULL,
    min_observed = 2,
    verbose = TRUE,
    ...
) {
    cluster_method <- match.arg(cluster_method)
    if (is.null(make_col)) make_col <- tti_default_make_col

    say <- function(...) if (isTRUE(verbose)) message(...)

    # With metadata the caller works in their own sample names throughout:
    # the internal column encoding is applied here and undone on the way out.
    if (!is.null(metadata)) {
        return(tti_run_from_metadata(
            dat, metadata, sample_col, subject_col, time_col, taxon_col,
            K, cluster_method, use_outliers, seed, min_observed, verbose
        ))
    }

    info <- tti_survey_input(
        dat, taxon_col, times, parse_fun, make_col, min_observed, say
    )
    dat_full <- tti_add_missing_columns(dat, info, say)

    res <- tti_run_pipeline(
        dat_full, taxon_col, parse_fun, info,
        K, cluster_method, use_outliers, seed, say
    )

    tti_run_result(res, info, taxon_col)
})


#' Print a tti_run object
#'
#' @param x An object of class \code{"tti_run"}.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.tti_run <- function(x, ...) {
    cat("TaxaTimeImpute run\n")
    cat("  taxa             :", nrow(x$completed), "\n")
    cat("  subjects         :", nrow(x$observed), "\n")
    cat("  missing samples  :", nrow(x$missing), "\n")
    cat("  cells imputed    :", nrow(x$imputed) - x$n_failed,
        "of", nrow(x$imputed), "\n")
    if (x$n_failed > 0) {
        cat("  cells left NA    :", x$n_failed, "\n")
    }
    thin <- sum(x$observed$n_observed < 2)
    if (thin > 0) {
        cat("  subjects with < 2 observations:", thin, "\n")
    }
    cat("\nUse $completed for the filled table, $imputed for the long form.\n")
    invisible(x)
}
