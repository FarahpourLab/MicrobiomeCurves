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
#' @param dat An abundance table: taxa in the \strong{row names}, one column
#'   per sample, named however the study names them. Values are either
#'   CLR-transformed abundances or raw counts, according to
#'   \code{abundance_type}.
#' @param metadata A data.frame with one row per sample, saying which subject
#'   each sample came from and when.
#' @param sample_col Character. Name of the \code{metadata} column holding
#'   sample identifiers, matching the columns of \code{dat}.
#' @param subject_col Character. Name of the \code{metadata} column
#'   identifying whom or what the repeated measurements were taken on: a
#'   mouse, a participant, a plot, a bioreactor. Nothing here assumes a
#'   species.
#' @param time_col Character. Name of the \code{metadata} column holding the
#'   time point. Numbers are used as they stand, so real spacing is
#'   respected. Labels such as \code{"baseline"} and \code{"week1"} are
#'   placed in order at equal spacing, taking the order from the factor's
#'   levels if it is a factor and from the order the rows appear otherwise.
#'   That assumption is reported, because the spacing changes every fitted
#'   curve.
#' @param abundance_type Either \code{"clr"}, meaning \code{dat} already
#'   holds centred log-ratios, or \code{"raw"}, meaning counts or relative
#'   abundances that are CLR-transformed here.
#' @param pseudocount Used only when \code{abundance_type = "raw"}. Either
#'   \code{"auto"}, replacing zeros per sample with a fraction of that
#'   sample's smallest non-zero value, or a single positive number added to
#'   every entry.
#' @param out_dir Optional directory to write results into. When given, the
#'   completed table is written as \code{imputed_abundance.tsv} and a record
#'   of the run as \code{imputation_log.txt}, holding the time points in
#'   order, the counts of samples, subjects and time points, every
#'   subject-timepoint that was missing and why, and any warning raised.
#' @param K Integer or \code{NULL}. Number of trajectory clusters, passed
#'   straight to \code{\link{tti_fit}}. \code{NULL} selects K per taxon by mean
#'   silhouette width.
#' @param cluster_method Character, passed to \code{\link{tti_fit}}.
#' @param use_outliers Logical, passed to \code{\link{tti_fit}}. \code{TRUE}
#'   enables functional-depth outlier screening.
#' @param seed Integer random seed, passed to \code{\link{tti_fit}}.
#' @param min_observed Integer. Observation count below which a subject is
#'   flagged in a warning. Default \code{2}.
#' @param verbose Logical. Print progress. Default \code{TRUE}.
#'
#' @return An object of class \code{"tti_run"}: a list with
#' \itemize{
#'   \item \code{completed}: the completed table, under the caller's own
#'     sample names, ordered by subject and then time. Observed values are
#'     unchanged. A column is added only for a subject-timepoint that had no
#'     sample; it is named from its subject and time.
#'   \item \code{metadata}: the metadata with a row for each created sample,
#'     carrying \code{imputed = TRUE}, so imputed samples stay
#'     distinguishable from measured ones.
#'   \item \code{design}: the \code{tti_design} the run was built from.
#'   \item \code{imputed}: long table of every imputed cell.
#'   \item \code{missing}: what was missing, and why.
#'   \item \code{observed}: observed time-point count per subject.
#'   \item \code{n_failed}: number of cells FPCA could not impute.
#'   \item \code{warnings}: every warning raised during the run.
#'   \item \code{files}: the paths written, when \code{out_dir} was given.
#'   \item \code{fit}: the underlying \code{"tti_fit"} object.
#' }
#'
#' @examples
#' # The bundled example in study form: taxa in the row names, samples named
#' # however the study names them, and a metadata sheet.
#' demo <- tti_demo_data()
#'
#' demo$counts[1:3, 1:4]
#' head(demo$metadata, 3)
#'
#' run <- suppressWarnings(tti_run(
#'     demo$counts, demo$metadata,
#'     sample_col = "sample", subject_col = "subject", time_col = "time",
#'     K = 1, verbose = FALSE
#' ))
#'
#' # what was missing, and why
#' run$missing
#'
#' # results under the original sample names
#' run$completed[1:3, 1:5]
#'
#' # samples created for subject-timepoints that had none
#' run$metadata[run$metadata$imputed, ]
#'
#' # Writing the completed table and a record of the run to a directory.
#' out <- file.path(tempdir(), "imputation")
#' run2 <- suppressWarnings(tti_run(
#'     demo$counts, demo$metadata,
#'     sample_col = "sample", subject_col = "subject", time_col = "time",
#'     out_dir = out, K = 1, verbose = FALSE
#' ))
#' basename(run2$files)
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
    metadata,
    sample_col,
    subject_col,
    time_col,
    abundance_type = c("clr", "raw"),
    pseudocount = "auto",
    out_dir = NULL,
    K = NULL,
    cluster_method = c("fpca", "kmeans_fd"),
    use_outliers = TRUE,
    seed = 123,
    min_observed = 2,
    verbose = TRUE,
    ...
) {
    tti_run_from_metadata(
        dat, metadata, sample_col, subject_col, time_col,
        match.arg(abundance_type), pseudocount, out_dir,
        K, match.arg(cluster_method), use_outliers, seed,
        min_observed, verbose
    )
})

#' @rdname tti_run
#' @export
setMethod("tti_run", "matrix", function(
    dat,
    metadata,
    sample_col,
    subject_col,
    time_col,
    abundance_type = c("clr", "raw"),
    pseudocount = "auto",
    out_dir = NULL,
    K = NULL,
    cluster_method = c("fpca", "kmeans_fd"),
    use_outliers = TRUE,
    seed = 123,
    min_observed = 2,
    verbose = TRUE,
    ...
) {
    # A matrix is the natural shape for taxa-in-row-names abundance data, so
    # it is accepted alongside a data.frame rather than being coerced by the
    # caller.
    tti_run_from_metadata(
        dat, metadata, sample_col, subject_col, time_col,
        match.arg(abundance_type), pseudocount, out_dir,
        K, match.arg(cluster_method), use_outliers, seed,
        min_observed, verbose
    )
})

#' Impute a table whose columns already encode subject and time
#'
#' @description
#' The layout the fitting code works on, with sample columns named
#' `"<subject>.<time>"`. Used by the `SummarizedExperiment` method and by the
#' metadata path once conversion has happened. Callers preparing their own
#' data use [tti_run()] with a metadata table instead.
#'
#' @param dat A data.frame in that layout.
#' @param taxon_col Character name of the taxon identifier column.
#' @param K,cluster_method,use_outliers,seed Passed to the fit.
#' @param times Optional numeric vector of the intended time points.
#' @param parse_fun,make_col Optional column-name parser and builder.
#' @param min_observed Integer. Subjects with fewer observed time points are
#'   reported.
#' @param verbose Logical. Whether to report progress.
#'
#' @return An object of class `tti_run`.
#'
#' @keywords internal
#' @noRd
tti_run_wide <- function(dat, taxon_col = "OTU_ID", K = NULL,
                         cluster_method = "fpca", use_outliers = TRUE,
                         seed = 123, times = NULL, parse_fun = NULL,
                         make_col = NULL, min_observed = 2, verbose = TRUE) {
    if (is.null(make_col)) make_col <- tti_default_make_col
    say <- function(...) if (isTRUE(verbose)) message(...)

    info <- tti_survey_input(
        dat, taxon_col, times, parse_fun, make_col, min_observed, say
    )
    dat_full <- tti_add_missing_columns(dat, info, say)

    res <- tti_run_pipeline(
        dat_full, taxon_col, parse_fun, info,
        K, cluster_method, use_outliers, seed, say
    )

    tti_run_result(res, info, taxon_col)
}


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
