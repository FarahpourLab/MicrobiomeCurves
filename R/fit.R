#' Flag outlying trajectories by functional depth
#'
#' @description
#' Fits FPCA to the supplied trajectories and marks those whose Fraiman-Muniz
#' depth falls below the `alpha` quantile.
#'
#' @param Ly List of numeric vectors, one per subject, holding the observed
#'   values.
#' @param Lt List of numeric vectors, one per subject, holding the observation
#'   times matching `Ly`.
#' @param alpha Numeric. Depth quantile below which a trajectory is flagged.
#'
#' @return Logical vector with one element per element of `Ly`. `TRUE` marks an
#'   outlier. Subjects that FPCA could not use are `FALSE`.
#'
#' @keywords internal
#' @noRd
detect_outliers_depth <- function(Ly, Lt, alpha = 0.05) {

    if (!requireNamespace("fda.usc", quietly = TRUE)) {
        stop("Package 'fda.usc' required for depth-based outliers")
    }

    fp <- mc_safe_fpca(Ly, Lt)
    if (is.null(fp)) return(rep(FALSE, length(Ly)))

    curves <- fitted(fp)
    grid <- fp$workGrid

    fdata_obj <- fda.usc::fdata(curves, argvals = grid)

    depth_obj <- fda.usc::depth.FM(fdata_obj)
    depth_vals <- depth_obj$dep

    cutoff <- quantile(depth_vals, alpha, na.rm = TRUE)

    outliers <- depth_vals < cutoff

    # mc_safe_fpca() fits only the subjects with at least two observations, so
    # `outliers` is indexed against those, not against every subject in Ly. Map
    # it back to full length; a subject FPCA could not use is not an outlier.
    #
    # When every subject was usable this is the identity, which is the case the
    # benchmark drivers always produce. Without it, the caller's
    # `names(outliers) <- names(Ly)` fails on any sparser input.
    obs_count <- vapply(Ly, function(x) sum(!is.na(x)), numeric(1))
    keep <- which(obs_count >= 2)

    if (length(outliers) == length(Ly)) {
        return(outliers)
    }

    out <- rep(FALSE, length(Ly))
    if (length(outliers) == length(keep)) {
        out[keep] <- outliers
    }

    return(out)
}

#' Silhouette-based selection of cluster number
#'
#' @param scores Numeric matrix of clustering features.
#' @param taxon_name Character taxon name used in the plot title.
#' @param max_K Maximum number of clusters to evaluate.
#' @param seed Random seed.
#' @param plot Logical. If TRUE, a silhouette diagnostic plot is attached to
#'   the result as the \code{"silhouette_plot"} attribute. It is attached
#'   rather than drawn, so that callers decide whether and where to render it.
#'
#' @return Integer optimal number of clusters, optionally carrying a
#'   \code{"silhouette_plot"} attribute.
#' @keywords internal
#' @noRd
select_K_silhouette_plot <- function(
    scores, taxon_name, max_K = 8, seed = 123, plot = TRUE) {
    if (!requireNamespace("cluster", quietly = TRUE)) {
        stop("Package 'cluster' required")
    }

    n <- nrow(scores)
    if (n < 3) return(1)

    max_K <- min(max_K, n - 1)
    if (max_K < 2) return(1)

    sil_scores <- rep(NA_real_, max_K)

    dmat <- stats::dist(scores)

    for (k in seq(2, max_K)) {
        mc_reset_rng(seed)
        km <- stats::kmeans(scores, centers = k, nstart = 10)
        sil <- cluster::silhouette(km$cluster, dmat)
        sil_scores[k] <- mean(sil[, 3])
    }

    valid_k <- which(!is.na(sil_scores))
    if (length(valid_k) == 0) return(1)

    K_opt <- which.max(sil_scores)

    if (plot) {
        df <- data.frame(k = seq_len(max_K), silhouette = sil_scores)

        p <- ggplot2::ggplot(df, ggplot2::aes(x = k, y = silhouette)) +
            ggplot2::geom_line(na.rm = TRUE) +
            ggplot2::geom_point(na.rm = TRUE) +
            ggplot2::geom_vline(
                xintercept = K_opt, linetype = "dashed", color = "red"
            ) +
            ggplot2::ggtitle(paste("Silhouette -", taxon_name)) +
            ggplot2::xlab("Number of clusters") +
            ggplot2::ylab("Mean silhouette width") +
            ggplot2::theme_minimal()

        attr(K_opt, "silhouette_plot") <- p
    }

    K_opt
}
#' @importFrom dplyr mutate filter left_join rename all_of
#' @importFrom tidyr expand_grid pivot_longer
#' @importFrom stats kmeans
NULL

#' Fit FPCA-based imputation model for longitudinal microbiome data
#'
#' @description
#' Fits a Functional Principal Component Analysis (FPCA) model to longitudinal
#' microbiome data and imputes missing whole-sample observations (i.e., all taxa
#' missing at a given subject-timepoint).
#'
#' @details
#' This function is the core modeling step of the MicrobiomeCurves workflow.
#'
#' For each taxon (feature), the function:
#' \enumerate{
#'   \item Masks specified subject-timepoint observations (entire samples)
#'   \item Constructs subject-specific trajectories
#'   \item Groups subjects using taxon-specific functional clustering (optional)
#'   \item Fits FPCA models within each cluster
#'   \item Predicts missing values using the fitted FPCA model
#' }
#'
#' Missingness is assumed to occur at the level of entire samples (i.e., all
#'  taxa
#' for a subject at a given time point are missing).
#'
#' @param prep A list returned by \code{\link{mc_prepare}} containing
#'   parsed data, replicate-time mapping, and masking information.
#' @param K Integer. Number of clusters for grouping subject trajectories.
#'   If \code{K = 1}, no clustering is performed and a single FPCA model is
#'   used.
#' @param cluster_method Character. How subjects are grouped once the taxon
#'   has been fitted. Both routes start from the same FPCA fit.
#'   \describe{
#'     \item{\code{"fpca"}}{The default. Runs k-means on the standardised
#'       FPC score vectors, comparing subjects through the leading modes of
#'       variation. This is the method behind the published results.}
#'     \item{\code{"kmeans_fd"}}{Runs functional k-means
#'       (\code{\link[fda.usc]{kmeans.fd}}) on the fitted curves, comparing
#'       whole trajectories under an L2 metric rather than a truncated score
#'       vector. \code{kmeans.fd()} fails on small sets of curves, which are
#'       common here once outliers are dropped; those cases fall back to
#'       \code{"fpca"} and are counted in the warning issued at the end of
#'       the fit.}
#'   }
#' @param use_outliers Logical. If \code{TRUE} (the default), trajectories are
#'   screened with Fraiman-Muniz functional depth and those below the 5\%
#'   quantile are dropped before the model is refitted. The target subject is
#'   always kept, even when flagged. If fewer than two trajectories would
#'   remain, screening is skipped for that taxon.
#' @param seed Integer. Random seed for reproducibility of clustering.
#'
#' @return
#' An object of class \code{"mc_fit"} containing:
#' \itemize{
#'   \item \code{dat}: Data with masked values applied
#'   \item \code{dat_orig}: Original unmodified data
#'   \item \code{taxon_col}: Taxon identifier column name
#'   \item \code{col_map}: Mapping of columns to replicate/time
#'   \item \code{reps}: Replicate identifiers
#'   \item \code{times}: Time points
#'   \item \code{clusters}: List of cluster assignments for each taxon
#'   \item \code{pred_long}: Data frame with true vs imputed values and errors
#'   \item \code{K}: Number of clusters used
#'   \item \code{seed}: Random seed used
#' }
#'
#' @examples
#' data(taxa_demo)
#'
#' prep <- mc_prepare(
#'     dat = taxa_demo,
#'     taxon_col = "OTU_ID",
#'     mask_list = data.frame(rep = "S01", time = 3)
#' )
#'
#' # K = 1 pools all subjects. K = NULL picks K per taxon by silhouette
#' # width, which is slower.
#' fit <- suppressWarnings(mc_fit(prep, K = 1))
#'
#' class(fit)
#' head(fit$pred_long[, c("species", "rep", "time", "imputed_value")])
#'
#' @seealso
#' \code{\link{mc_prepare}},
#' \code{\link{mc_impute}},
#' \code{\link{mc_run}}
#'
#' @export
mc_fit <- function(prep, K = NULL,
                    cluster_method = c("fpca", "kmeans_fd"),
                    use_outliers = TRUE, seed = 123) {
    # Reproducible for a given seed, but the caller's own random stream is
    # restored when this function returns. See R/rng.R.
    mc_preserve_rng()
    mc_reset_rng(seed)

    mc_check_prep(prep)

    dat_orig <- prep$dat
    taxon_col <- prep$taxon_col
    reps <- prep$reps
    cluster_method <- match.arg(cluster_method)
    species_vec <- dat_orig[[taxon_col]]

    pred_long <- mc_build_targets(
        dat_orig, taxon_col, prep$mask_pairs, nrow(dat_orig)
    )

    reshaped <- mc_mask_and_reshape(dat_orig, prep$col_map, prep$mask_pairs)
    dat <- reshaped$dat
    long <- reshaped$long

    # Count FPCA outcomes across the whole run so they can be reported once
    # rather than at each of the many thousands of individual fits.
    mc_diag_reset()
    on.exit(mc_diag_clear(), add = TRUE)

    ft <- mc_fit_all_taxa(
        long, taxon_col, species_vec, reps, pred_long,
        use_outliers, K, seed, cluster_method
    )
    pred_long <- dplyr::mutate(
        ft$pred_long, se = (true_value - imputed_value)^2
    )

    mc_diag_report(nrow(pred_long))
    mc_warn_unimputed(pred_long)

    mc_fit_result(
        dat, dat_orig, prep, ft, pred_long, K, cluster_method,
        use_outliers, seed
    )
}

#' Assemble the object returned by mc_fit()
#'
#' @param dat The table with masked cells set to `NA`.
#' @param dat_orig The table as supplied.
#' @param prep The object from [mc_prepare()].
#' @param ft List returned by `mc_fit_all_taxa()`.
#' @param pred_long Prediction table carrying the squared errors.
#' @param K Integer or `NULL`, the requested number of clusters.
#' @param cluster_method Character, the clustering route that was used.
#' @param use_outliers Logical, whether trajectories were screened.
#' @param seed Integer random seed.
#'
#' @return An object of class `mc_fit`.
#'
#' @keywords internal
#' @noRd
mc_fit_result <- function(dat, dat_orig, prep, ft, pred_long, K,
                           cluster_method, use_outliers, seed) {
    structure(
        c(
            list(dat = dat, dat_orig = dat_orig),
            prep[c("taxon_col", "col_map", "reps", "times")],
            list(
                clusters = ft$clusters, pred_long = pred_long,
                K = K, cluster_method = cluster_method,
                outliers = ft$outliers, use_outliers = use_outliers,
                seed = seed
            )
        ),
        class = "mc_fit"
    )
}
