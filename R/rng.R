#' Set the RNG stream without leaving the caller's stream disturbed
#'
#' @description
#' Internal replacement for a direct `set.seed()` call.
#'
#' @details
#' The clustering step draws random starting centres, so a `seed` argument is
#' part of the user-facing API and results must be reproducible for a given
#' seed. Calling `set.seed()` inside a package function would achieve that but
#' would also silently reset the caller's own random stream, which is why
#' Bioconductor asks packages not to do it.
#'
#' `mc_reset_rng()` produces exactly the RNG state that `set.seed(seed)`
#' produces, including the stream left behind for subsequent draws, and
#' `mc_preserve_rng()` restores whatever stream the caller had once the
#' calling function returns. Together they keep results identical to earlier
#' versions of this code while leaving the caller's stream untouched across
#' the call.
#'
#' @param seed Integer seed.
#'
#' @return `NULL`, invisibly. Called for its effect on the RNG state.
#'
#' @keywords internal
#' @noRd
mc_reset_rng <- function(seed) {
    state <- withr::with_seed(
        seed,
        get(".Random.seed", envir = globalenv())
    )
    assign(".Random.seed", state, envir = globalenv())
    invisible(NULL)
}

#' Load the bundled demo table
#'
#' @description
#' `LazyData` is switched off, as Bioconductor asks, so `taxa_demo` is not
#' attached automatically. Internal code loads it through here rather than
#' relying on lazy-loading.
#'
#' @return The `taxa_demo` data frame.
#'
#' @keywords internal
#' @noRd
mc_demo_table <- function() {
    env <- new.env(parent = emptyenv())
    utils::data("taxa_demo", package = "MicrobiomeCurves", envir = env)
    env$taxa_demo
}

#' Restore the caller's RNG stream when the calling function exits
#'
#' @param frame Environment whose exit triggers the restore.
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
#' @noRd
mc_preserve_rng <- function(frame = parent.frame()) {
    withr::local_preserve_seed(.local_envir = frame)
    invisible(NULL)
}
