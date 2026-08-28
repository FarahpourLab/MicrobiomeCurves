############################################
# Analytic CI (internal)
############################################

#' Analytic confidence interval for FPCA prediction
#'
#' @description
#' Computes analytic (model-based) confidence intervals for FPCA predictions
#' at a given time point for a single subject trajectory.
#'
#' @details
#' This function uses the FPCA model outputs (mean function, eigenfunctions,
#' eigenvalues, and noise variance) to compute the conditional distribution
#' of the predicted trajectory.
#'
#' The prediction is based on the Best Linear Unbiased Predictor (BLUP)
#' of the subject-specific scores, and uncertainty is derived from the
#' posterior variance of the FPCA scores.
#'
#' The variance of the predicted value is given by:
#' \deqn{
#' Var(\hat{Y}(t)) = \phi(t)^T V_{\xi} \phi(t) + \sigma^2
#' }
#'
#' where \eqn{V_{\xi}} is the conditional covariance of the FPCA scores.
#'
#' @param fp An FPCA object returned by \code{fdapace::FPCA}.
#'
#' @param obs_times Numeric vector of observed time points for the target
#'  subject.
#'
#' @param obs_values Numeric vector of observed values corresponding to
#'  \code{obs_times}.
#'
#' @param pred_time Numeric scalar. Time point at which prediction is computed.
#'
#' @param include_noise Logical. If \code{TRUE}, includes measurement noise
#'   variance in the prediction variance. Default is \code{TRUE}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{mean}: Predicted value at \code{pred_time}
#'   \item \code{se}: Standard error of the prediction
#'   \item \code{lower}: Lower bound of 95\% confidence interval
#'   \item \code{upper}: Upper bound of 95\% confidence interval
#' }
#'
#' @seealso
#' \code{\link{mc_ci}}
#'
#' @keywords internal
#' @noRd
mc_analytic_ci <- function(fp, obs_times, obs_values, pred_time,
                            include_noise = TRUE) {
    if (is.null(fp)) return(NULL)

    parts <- mc_fpca_parts(fp)
    if (is.null(parts)) return(NULL)

    grid <- parts$grid
    mu <- parts$mu
    phi <- parts$phi
    lambda <- parts$lambda
    sigma2 <- parts$sigma2

    # vapply, not sapply: for an empty obs_times sapply() returns an empty list,
    # and `phi[list(), , drop = FALSE]` fails with "invalid subscript type".
    # A subject can legitimately have no observations left once its samples are
    # masked or missing. vapply() yields integer(0) there, and the algebra below
    # then reduces to the population mean with the prior score variance, which
    # is
    # the correct answer when nothing subject-specific is known.
    # For non-empty obs_times the two are identical.
    obs_idx <- vapply(
        obs_times, function(tt) which.min(abs(grid - tt)), integer(1)
    )
    pred_idx <- which.min(abs(grid - pred_time))

    blup <- mc_score_blup(
        Phi_i = phi[obs_idx, , drop = FALSE],
        y_center = obs_values - mu[obs_idx],
        lambda = lambda,
        sigma2 = sigma2
    )
    if (is.null(blup)) return(NULL)

    phi_t <- phi[pred_idx, ]
    pred_mean <- as.numeric(mu[pred_idx] + crossprod(phi_t, blup$xi_hat))

    pred_var <- as.numeric(t(phi_t) %*% blup$V_xi %*% phi_t)
    if (include_noise) pred_var <- pred_var + sigma2
    pred_var <- max(pred_var, 0)
    pred_sd <- sqrt(pred_var)

    list(
        mean = pred_mean,
        se = pred_sd,
        lower = pred_mean - 1.96 * pred_sd,
        upper = pred_mean + 1.96 * pred_sd
    )
}

#' Best linear unbiased predictor of the FPCA scores
#'
#' @param Phi_i Matrix of eigenfunctions at the observed times, one row per
#'   observation.
#' @param y_center Numeric vector of observed values minus the mean function.
#' @param lambda Numeric vector of eigenvalues.
#' @param sigma2 Numeric measurement-error variance.
#'
#' @return List with the score covariance `V_xi` and the predicted scores
#'   `xi_hat`, or `NULL` when the system cannot be solved.
#'
#' @keywords internal
#' @noRd
mc_score_blup <- function(Phi_i, y_center, lambda, sigma2) {
    Lambda_inv <- diag(1 / lambda, nrow = length(lambda))

    V_xi <- tryCatch(
        solve(Lambda_inv + crossprod(Phi_i) / sigma2),
        error = function(e) NULL
    )
    if (is.null(V_xi)) return(NULL)

    list(
        V_xi = V_xi,
        xi_hat = V_xi %*% crossprod(Phi_i, y_center) / sigma2
    )
}


############################################
# Bootstrap CI (internal)
############################################


#' Bootstrap confidence intervals for FPCA trajectories
#'
#' @description
#' Computes bootstrap-based confidence intervals for FPCA-fitted trajectories
#' in longitudinal microbiome data by resampling subjects (replicates).
#'
#' @details
#' This function performs nonparametric bootstrap by resampling subjects with
#' replacement and refitting the FPCA model at each iteration.
#'
#' For each bootstrap sample:
#' \enumerate{
#'   \item Subjects are resampled with replacement
#'   \item The FPCA model is refit using the resampled data
#'   \item The trajectory for the target subject is extracted
#' }
#'
#' Pointwise confidence intervals are constructed using empirical quantiles
#' (2.5\% and 97.5\%) across bootstrap samples.
#'
#' The target subject is always included in each bootstrap iteration to ensure
#' stable prediction.
#'
#' @param Ly A list of numeric vectors containing observed values for each
#'  subject.
#'
#' @param Lt A list of numeric vectors containing observation times
#'  corresponding
#'   to \code{Ly}.
#'
#' @param target_rep Character string specifying the subject (replicate) for
#'  which
#'   confidence intervals are computed.
#'
#' @param B Integer. Number of bootstrap samples. Default is \code{500}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{grid}: Evaluation time grid used in FPCA
#'   \item \code{lower}: Lower confidence bound (2.5\% quantile)
#'   \item \code{upper}: Upper confidence bound (97.5\% quantile)
#' }
#'
#' @seealso
#' \code{\link{mc_ci}}
#'
#' @keywords internal
#' @noRd
mc_bootstrap_ci <- function(Ly, Lt, target_rep, B = 500) {
    fp0 <- tryCatch(
        fdapace::FPCA(Ly = Ly, Lt = Lt, optns = list(dataType = "Sparse")),
        error = function(e) NULL
    )
    if (is.null(fp0)) return(NULL)

    grid <- fp0$workGrid
    members <- names(Ly)
    curves <- matrix(NA, nrow = B, ncol = length(grid))

    for (b in seq_len(B)) {
        samp <- sample(members, length(members), replace = TRUE)

        Ly_b <- Ly[samp]
        Lt_b <- Lt[samp]

        names(Ly_b) <- paste0(samp, "_", b)
        names(Lt_b) <- names(Ly_b)

        Ly_b[[target_rep]] <- Ly[[target_rep]]
        Lt_b[[target_rep]] <- Lt[[target_rep]]

        fp_b <- tryCatch(
            fdapace::FPCA(
                Ly = Ly_b, Lt = Lt_b,
                optns = list(dataType = "Sparse")
            ),
            error = function(e) NULL
        )

        if (!is.null(fp_b)) {
            fit_b <- fitted(fp_b)
            idx <- match(target_rep, names(Ly_b))
            if (!is.na(idx)) curves[b, ] <- fit_b[idx, ]
        }
    }

    lower <- apply(curves, 2, stats::quantile, 0.025, na.rm = TRUE)
    upper <- apply(curves, 2, stats::quantile, 0.975, na.rm = TRUE)

    list(grid = grid, lower = lower, upper = upper)
}

############################################
# Main user function
############################################


#' Compute uncertainty bands for FPCA imputation
#'
#' @description
#' Computes confidence intervals for FPCA-based imputed trajectories in
#' longitudinal microbiome data using analytic, bootstrap, or combined methods.
#'
#' @details
#' This function estimates uncertainty for imputed values at a given
#' subject-timepoint by fitting a cluster-specific FPCA model and deriving
#' confidence bands along the time grid.
#'
#' Two types of uncertainty estimation are supported:
#' \itemize{
#'   \item \strong{Analytic}: Model-based confidence intervals derived from
#'   the conditional variance of FPCA scores (fast, parametric)
#'
#'   \item \strong{Bootstrap}: Nonparametric confidence intervals obtained by
#'   resampling subjects and refitting FPCA (robust, data-driven)
#' }
#'
#' When \code{method = "both"}, both uncertainty estimates are returned.
#'
#' The FPCA model is fitted using only subjects within the same cluster as the
#' target subject, allowing borrowing of information from similar trajectories.
#'
#' @param fit An object of class \code{"mc_fit"} returned by
#'   \code{\link{mc_fit}}.
#'
#' @param species_name Character. Name of the taxon (feature) to analyze.
#'
#' @param rep_id Character. Target subject (replicate) for which uncertainty
#'   is computed.
#'
#' @param time_id Numeric. Time point at which the value is considered missing
#'   and imputed.
#'
#' @param method Character. Type of confidence interval to compute:
#'   \code{"analytic"}, \code{"bootstrap"}, or \code{"both"}.
#'
#' @param B Integer. Number of bootstrap samples used when
#'   \code{method = "bootstrap"} or \code{"both"}. Default is \code{200}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{fp}: FPCA object fitted on cluster members
#'   \item \code{analytic}: Data frame with analytic confidence bands
#'   (columns: \code{time}, \code{lower}, \code{upper}), or \code{NULL}
#'   \item \code{bootstrap}: Data frame with bootstrap confidence bands,
#'   or \code{NULL}
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
#' fit <- suppressWarnings(mc_fit(prep, K = 1))
#'
#' # analytic intervals are cheap
#' ci <- suppressWarnings(
#'     mc_ci(
#'         fit = fit,
#'         species_name = fit$pred_long$species[1],
#'         rep_id = "S01",
#'         time_id = 3,
#'         method = "analytic"
#'     )
#' )
#'
#' head(ci$analytic)
#'
#' # Bootstrap intervals refit the model B times. B is small here to keep the
#' # example quick; use a few hundred for real work.
#' ci_both <- suppressWarnings(
#'     mc_ci(
#'         fit = fit,
#'         species_name = fit$pred_long$species[1],
#'         rep_id = "S01",
#'         time_id = 3,
#'         method = "both",
#'         B = 25
#'     )
#' )
#'
#' head(ci_both$bootstrap)
#'
#' @seealso
#' \code{\link{mc_fit}},
#' \code{\link{mc_plot}}
#'
#' @export
mc_ci <- function(
    fit,
    species_name,
    rep_id,
    time_id,
    method = c("analytic", "bootstrap", "both"),
    B = 200
) {
    method <- match.arg(method)

    clusters <- fit$clusters[[species_name]]
    if (is.null(clusters)) {
        stop("No clustering found for this species")
    }

    cl_id <- clusters[rep_id]
    members <- names(clusters)[clusters == cl_id]

    W_sp <- mc_subject_time_matrix(fit, species_name)
    fpca <- mc_cluster_fpca(W_sp, members, rep_id, time_id, fit$times)
    fp <- fpca$fp

    out <- list(fp = fp, analytic = NULL, bootstrap = NULL)

    if (method %in% c("analytic", "both")) {
        out$analytic <- mc_analytic_band(
            fp,
            obs_times = fpca$Lt[[rep_id]],
            obs_values = fpca$Ly[[rep_id]]
        )
    }

    if (method %in% c("bootstrap", "both")) {
        band <- mc_bootstrap_ci(
            fpca$Ly, fpca$Lt,
            target_rep = rep_id, B = B
        )
        if (!is.null(band)) {
            out$bootstrap <- tibble::tibble(
                time = band$grid,
                lower = band$lower,
                upper = band$upper
            )
        }
    }

    out
}
