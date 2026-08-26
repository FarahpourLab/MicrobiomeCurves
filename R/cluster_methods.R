# The two ways of grouping subjects before a cell is imputed.
#
# Both start from the same FPCA fit of the taxon. They differ in what is
# handed to k-means:
#
#   "fpca"      clusters the standardised FPC score vectors, so subjects are
#               compared through the leading modes of variation. This is the
#               package default and the method behind the published results.
#
#   "kmeans_fd" clusters the fitted curves themselves with
#               fda.usc::kmeans.fd(), so subjects are compared by their whole
#               trajectory under an L2 metric rather than by a truncated
#               score vector.
#
# kmeans.fd() needs a reasonable number of curves and fails outright on very
# small sets, which is common here once outliers are dropped. Rather than
# lose the taxon, that case falls back to the score-based route and is
# counted so tti_fit() can report how often it happened.

#' Assign subjects to clusters
#'
#' @param fp Fitted FPCA object for the subjects in play.
#' @param scores Matrix of standardised FPC scores, one row per subject.
#' @param k_use Integer number of clusters to form.
#' @param subjects Character vector of subject identifiers, in the row order
#'   of `scores`.
#' @param method Either `"fpca"` or `"kmeans_fd"`.
#'
#' @return Named integer vector giving a cluster for each subject.
#'
#' @keywords internal
#' @noRd
tti_assign_clusters <- function(fp, scores, k_use, subjects, method) {
    if (identical(method, "kmeans_fd")) {
        cl <- tti_kmeans_fd(fp, k_use, subjects)
        if (!is.null(cl)) {
            return(cl)
        }
        tti_diag_bump("kfd_fallback")
    }

    km <- kmeans(scores, centers = k_use, nstart = 20)
    clusters <- km$cluster
    names(clusters) <- subjects
    clusters
}

#' Cluster the fitted curves with functional k-means
#'
#' @description
#' Builds an `fdata` object from the FPCA fitted trajectories on the shared
#' work grid and passes it to [fda.usc::kmeans.fd()]. Returns `NULL` when the
#' fit is unusable or `kmeans.fd()` fails, which it does on small sets of
#' curves, so the caller can fall back.
#'
#' @param fp Fitted FPCA object.
#' @param k_use Integer number of clusters to form.
#' @param subjects Character vector of subject identifiers, in the row order
#'   of the fitted curves.
#'
#' @return Named integer vector of cluster assignments, or `NULL`.
#'
#' @keywords internal
#' @noRd
tti_kmeans_fd <- function(fp, k_use, subjects) {
    if (!requireNamespace("fda.usc", quietly = TRUE)) {
        stop(
            "Package 'fda.usc' is required for cluster_method = ",
            "\"kmeans_fd\".",
            call. = FALSE
        )
    }

    curves <- tryCatch(stats::fitted(fp), error = function(e) NULL)
    if (is.null(curves) || !is.matrix(curves) ||
        nrow(curves) != length(subjects)) {
        return(NULL)
    }

    grid <- fp$workGrid
    if (is.null(grid) || ncol(curves) != length(grid)) {
        return(NULL)
    }

    rownames(curves) <- subjects

    res <- tryCatch(
        suppressWarnings(suppressMessages(
            fda.usc::kmeans.fd(
                fda.usc::fdata(curves, argvals = grid),
                ncl = k_use,
                draw = FALSE
            )
        )),
        error = function(e) NULL
    )

    if (is.null(res) || is.null(res$cluster)) {
        return(NULL)
    }

    clusters <- as.integer(res$cluster)
    if (length(clusters) != length(subjects) || anyNA(clusters)) {
        return(NULL)
    }

    names(clusters) <- subjects
    clusters
}
