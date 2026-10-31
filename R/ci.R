############################################
# Analytic CI (internal)
############################################

#' Prediction and interval for one subject at one time, from fdapace
#'
#' @description
#' The value is the one fdapace itself computed: its mean curve plus that
#' subject's own scores, read at the grid point nearest the target time. That
#' is exactly what `fitted(fp)` returns, so it matches the original analysis
#' scripts.
#'
#' The scores used to be recomputed here, from nearest grid points and a
#' floored noise variance. On real data that diverged: for a subject with as
#' many observed points as there were components and a noise variance of
#' 1e-6, the recomputed scores came out three orders of magnitude too large,
#' and the prediction with them. fdapace already returns the scores it used,
#' so they are taken from there.
#'
#' The interval comes from the same place, fdapace's per-subject score
#' covariance, so it is consistent with the value.
#'
#' @param fp Fitted FPCA object carrying the `mc_subjects` attribute.
#' @param subject Character identifier of the subject to predict for.
#' @param pred_time Numeric time point to predict at.
#' @param include_noise Add the noise variance to the prediction variance.
#'
#' @return List with `mean`, `se`, `lower` and `upper`, or `NULL` when the
#'   subject was not part of the fit. The interval is `NA` when fdapace
#'   returned no usable score covariance; the value is still given, because
#'   that is what the imputation needs.
#'
#' @keywords internal
#' @noRd
mc_analytic_ci <- function(fp, subject, pred_time, include_noise = TRUE) {
    if (is.null(fp)) {
        return(NULL)
    }

    grid <- fp$workGrid
    mu <- fp$mu
    phi <- fp$phi
    xi <- fp$xiEst
    if (is.null(grid) || is.null(mu) || is.null(phi) || is.null(xi)) {
        return(NULL)
    }

    subjects <- attr(fp, "mc_subjects")
    i <- if (is.null(subjects)) NA_integer_ else match(subject, subjects)

    idx <- which.min(abs(grid - pred_time))
    phi_t <- phi[idx, ]
    s2 <- fp$sigma2
    if (is.null(s2) || !is.finite(s2)) s2 <- 0

    if (is.na(i) || i > nrow(xi)) {
        # The subject is not in the fit, because it had too few observations
        # for fdapace to score it. Nothing subject-specific is known, so the
        # prediction is the population mean curve, which is the limit of the
        # score estimate when there is nothing to condition on. The spread is
        # the prior one, from the eigenvalues.
        #
        # The original scripts index fitted(fp) by the subject's position
        # among all subjects, which in this situation points at a different
        # subject's curve or past the end of the matrix. That is not
        # reproduced here.
        lambda <- fp$lambda
        prior <- if (is.null(lambda)) NA_real_ else sum(lambda * phi_t^2)
        v <- if (is.na(prior)) NA_real_ else prior + if (include_noise) s2 else 0
        sd0 <- if (is.na(v)) NA_real_ else sqrt(max(v, 0))
        return(list(
            mean = as.numeric(mu[idx]),
            se = sd0,
            lower = as.numeric(mu[idx]) - 1.96 * sd0,
            upper = as.numeric(mu[idx]) + 1.96 * sd0
        ))
    }

    # identical to fitted(fp)[i, idx]
    pred_mean <- as.numeric(mu[idx] + crossprod(phi_t, xi[i, ]))

    pred_sd <- NA_real_
    xv <- fp$xiVar
    if (!is.null(xv) && length(xv) >= i) {
        V <- xv[[i]]
        if (is.matrix(V) && all(dim(V) == length(phi_t)) &&
            all(is.finite(V))) {
            v <- as.numeric(t(phi_t) %*% V %*% phi_t)
            if (include_noise) {
                v <- v + s2
            }
            pred_sd <- sqrt(max(v, 0))
        }
    }

    list(
        mean = pred_mean,
        se = pred_sd,
        lower = pred_mean - 1.96 * pred_sd,
        upper = pred_mean + 1.96 * pred_sd
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
#' fit <- suppressWarnings(mc_fit(prep, C = 1))
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
        out$analytic <- mc_analytic_band(fp, subject = rep_id)
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

#' Nearest positive semi-definite version of a covariance matrix
#'
#' @description
#' A score covariance returned by fdapace is a posterior covariance and so is
#' positive semi-definite in exact arithmetic, but rounding can leave a tiny
#' negative eigenvalue, and `chol()` then fails. Clamping the eigenvalues at
#' zero gives the nearest matrix that factorises.
#'
#' @param A Square numeric matrix.
#'
#' @return `A`, symmetrised, with negative eigenvalues clamped.
#'
#' @keywords internal
#' @noRd
mc_psd <- function(A) {
    e <- eigen((A + t(A)) / 2, symmetric = TRUE)
    e$vectors %*% diag(pmax(e$values, 1e-10), length(e$values)) %*%
        t(e$vectors)
}

#' Interval bounds for one cell, read off simulated draws
#'
#' @description
#' The analytic interval is `mean +- 1.96 se`, which is symmetric and so can
#' reach below the detection floor of the sample it belongs to. This function
#' simulates the same posterior instead of solving it, holds every draw at or
#' above that floor, and reads the 2.5% and 97.5% points off the result. The
#' bounds that come back are asymmetric, by as much as the floor is binding.
#'
#' @details
#' The draws reproduce the analytic interval when the floor does not bite, so
#' nothing here is a different model: `mu + phi' xi` with
#' `xi ~ N(xiEst, xiVar)` and measurement noise `sigma2` is exactly the
#' distribution whose standard deviation [mc_analytic_ci()] computes in closed
#' form.
#'
#' The draws run on their own RNG stream. The fit sets one seed and the
#' clustering step consumes from it for the rest of the run, so taking random
#' numbers here without restoring the stream would change which subjects enter
#' later fits, and with them the imputed values.
#'
#' @param fp A fitted `fdapace::FPCA` object.
#' @param subject Character identifier of the target subject.
#' @param pred_time Numeric time point to predict.
#' @param floor_v Numeric detection floor of the cell's own sample. The bounds
#'   are `NULL` when this is not finite.
#' @param ndraw Integer number of draws.
#' @param cell_seed Integer seed for this cell, so the bounds are reproducible.
#'
#' @return List with `lower`, `upper`, and `share`, the fraction of draws the
#'   floor held up. `NULL` when the cell cannot be drawn from.
#'
#' @keywords internal
#' @noRd
mc_cell_bounds <- function(fp, subject, pred_time, floor_v,
                            ndraw = 10000L, cell_seed = 1L) {
    if (!is.finite(floor_v) || ndraw < 1L) {
        return(NULL)
    }

    grid <- fp$workGrid
    mu <- fp$mu
    phi <- fp$phi
    xi <- fp$xiEst
    if (is.null(grid) || is.null(mu) || is.null(phi) || is.null(xi)) {
        return(NULL)
    }

    idx <- which.min(abs(grid - pred_time))
    phi_t <- phi[idx, ]
    k <- length(phi_t)
    s2 <- fp$sigma2
    if (is.null(s2) || !is.finite(s2) || s2 < 0) s2 <- 0

    subjects <- attr(fp, "mc_subjects")
    i <- if (is.null(subjects)) NA_integer_ else match(subject, subjects)

    if (is.na(i) || i > nrow(xi)) {
        # Not scored by fdapace, so the draws come from the population, which
        # is the same fallback [mc_analytic_ci()] takes.
        lambda <- fp$lambda
        if (is.null(lambda) || length(lambda) != k) {
            return(NULL)
        }
        centre <- as.numeric(mu[idx])
        v_mat <- diag(lambda, k)
    } else {
        centre <- as.numeric(mu[idx] + crossprod(phi_t, xi[i, ]))
        xv <- fp$xiVar
        if (is.null(xv) || length(xv) < i) {
            return(NULL)
        }
        v_mat <- xv[[i]]
        if (!is.matrix(v_mat) || !all(dim(v_mat) == k) ||
            !all(is.finite(v_mat))) {
            return(NULL)
        }
    }

    fac <- tryCatch(t(chol(mc_psd(v_mat))), error = function(e) NULL)
    if (is.null(fac)) {
        return(NULL)
    }

    draws <- withr::with_preserve_seed({
        mc_reset_rng(cell_seed)
        z <- matrix(stats::rnorm(k * ndraw), k, ndraw)
        as.numeric(centre + crossprod(phi_t, fac %*% z)) +
            stats::rnorm(ndraw, 0, sqrt(s2))
    })

    share <- mean(draws < floor_v)
    if (share == 0) {
        # Nothing was censored, so the draws are estimating a distribution
        # that is symmetric, and their quantiles differ from the closed form
        # only by simulation error. Reporting that error would show as
        # asymmetry that is not there, so the exact bounds are used instead.
        sd_a <- sqrt(max(
            as.numeric(t(phi_t) %*% v_mat %*% phi_t) + s2, 0
        ))
        return(list(
            lower = centre - 1.96 * sd_a,
            upper = centre + 1.96 * sd_a,
            share = 0
        ))
    }

    draws <- pmax(draws, floor_v)
    q <- stats::quantile(draws, c(0.025, 0.975), names = FALSE)
    list(lower = q[1], upper = q[2], share = share)
}
