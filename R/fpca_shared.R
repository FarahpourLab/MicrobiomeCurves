# Helpers shared by mc_ci() and mc_plot().

#' Build the subject-by-time matrix for one taxon
#'
#' @param fit Object of class `mc_fit`.
#' @param species_name Character name of the taxon.
#'
#' @return Numeric matrix with subjects as rows and time points as columns.
#'
#' @keywords internal
#' @noRd
mc_subject_time_matrix <- function(fit, species_name) {
    dat <- fit$dat
    col_map <- fit$col_map
    reps <- fit$reps
    times <- fit$times
    taxon_col <- fit$taxon_col

    y_all <- as.numeric(dat[dat[[taxon_col]] == species_name, col_map$col])

    W_sp <- matrix(
        NA,
        nrow = length(reps), ncol = length(times),
        dimnames = list(reps, times)
    )
    for (k in seq_len(nrow(col_map))) {
        W_sp[col_map$rep[k], as.character(col_map$time[k])] <- y_all[k]
    }
    W_sp
}

#' Fit FPCA to the target subject's cluster, with the target point masked
#'
#' @param W_sp Subject-by-time matrix from [mc_subject_time_matrix()].
#' @param members Character vector of subjects in the target's cluster.
#' @param rep_id Character identifier of the target subject.
#' @param time_id Numeric time point being imputed.
#' @param times Numeric vector of all time points.
#'
#' @return List with the fitted FPCA object `fp`, its `grid`, the fitted
#'   curve matrix `fit_mat`, and the `Ly` and `Lt` lists used for the fit.
#'
#' @keywords internal
#' @noRd
mc_cluster_fpca <- function(W_sp, members, rep_id, time_id, times) {
    Lt <- Ly <- vector("list", length(members))
    names(Lt) <- names(Ly) <- members

    for (r in members) {
        yy <- W_sp[r, ]
        if (r == rep_id) {
            yy[which(times == time_id)] <- NA
        }
        keep <- which(!is.na(yy))
        Lt[[r]] <- times[keep]
        Ly[[r]] <- yy[keep]
    }

    fp <- fdapace::FPCA(
        Ly = Ly, Lt = Lt,
        optns = list(dataType = "Sparse")
    )

    list(
        fp = fp, grid = fp$workGrid, fit_mat = fitted(fp),
        Ly = Ly, Lt = Lt
    )
}


#' Analytic confidence band along the whole time grid
#'
#' @param fp Fitted FPCA object.
#' @param obs_times Numeric vector of observed times for the target subject.
#' @param obs_values Numeric vector of observed values for the target subject.
#'
#' @return Tibble with `time`, `lower` and `upper`, or `NULL` when no point on
#'   the grid could be computed.
#'
#' @keywords internal
#' @noRd
mc_analytic_band <- function(fp, obs_times, obs_values) {
    grid <- fp$workGrid

    ci_list <- lapply(grid, function(tt) {
        mc_analytic_ci(fp, obs_times, obs_values, tt, include_noise = TRUE)
    })

    valid <- !vapply(ci_list, is.null, logical(1))
    if (!any(valid)) {
        return(NULL)
    }

    tibble::tibble(
        time = grid[valid],
        lower = vapply(ci_list[valid], function(x) x$lower, numeric(1)),
        upper = vapply(ci_list[valid], function(x) x$upper, numeric(1))
    )
}

#' Pull the pieces of a fitted FPCA model needed for prediction
#'
#' @param fp An FPCA object from `fdapace::FPCA`.
#'
#' @return List with `grid`, `mu`, `phi`, `lambda` and `sigma2`, with
#'   the variances floored away from zero, or `NULL` when the model is
#'   unusable.
#'
#' @keywords internal
#' @noRd
mc_fpca_parts <- function(fp) {
    grid <- fp$workGrid
    mu <- fp$mu
    phi <- fp$phi
    lambda <- fp$lambda
    sigma2 <- fp$sigma2

    if (is.null(grid) || is.null(mu) || is.null(phi) ||
        is.null(lambda)) {
        return(NULL)
    }
    if (length(lambda) == 0) {
        return(NULL)
    }

    list(
        grid = grid, mu = mu, phi = phi,
        lambda = pmax(lambda, 1e-8),
        sigma2 = ifelse(
            is.null(sigma2) || is.na(sigma2), 1e-6, max(sigma2, 1e-6)
        )
    )
}
