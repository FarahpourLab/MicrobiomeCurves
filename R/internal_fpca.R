#' Safe wrapper for FPCA fitting
#'
#' @description
#' Internal utility function that safely fits a Functional Principal Component
#'  Analysis (FPCA)
#' model using the \pkg{fdapace} package. It filters out subjects with
#'  insufficient observations
#' and returns \code{NULL} if model fitting is not feasible.
#'
#' @details
#' The function requires at least two subjects with at least two observed time
#'  points each.
#' If these conditions are not met, or if the FPCA model fails to converge, the
#'  function
#' returns \code{NULL} instead of throwing an error.
#'
#' @param Ly A list of numeric vectors containing observed values for each
#'  subject.
#'
#' @param Lt A list of numeric vectors containing observation times
#'  corresponding to \code{Ly}.
#'
#' @return
#' An FPCA object returned by \code{fdapace::FPCA}, or \code{NULL} if fitting
#'  fails.
#'
#' @keywords internal
#'
#' @noRd
tti_safe_fpca <- function(Ly, Lt) {
    obs_count <- vapply(Ly, function(x) sum(!is.na(x)), numeric(1))
    keep <- which(obs_count >= 2)

    tti_diag_bump("calls")

    if (length(keep) < 2) {
        tti_diag_bump("too_few")
        return(NULL)
    }

    Ly_sub <- Ly[keep]
    Lt_sub <- Lt[keep]

    if (length(keep) <= 3) {
        tti_diag_bump("small")
    }

    # fdapace notes small samples and reset options on every call. Those
    # notes are counted above and reported once by tti_fit(), so they are
    # muffled here rather than repeated for each of the many thousands of
    # fits a real table produces.
    fit <- suppressMessages(tryCatch(
        fdapace::FPCA(
            Ly = Ly_sub,
            Lt = Lt_sub,
            optns = list(dataType = "Sparse")
        ),
        error = function(e) NULL
    ))

    if (is.null(fit)) {
        tti_diag_bump("failed")
    }
    fit
}
