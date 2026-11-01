# Drawing per-taxon uncertainty.
#
# One image per taxon, each imputed value shown with its interval. An image
# is readable on its own, and the set as a whole is the thing to look at
# before treating imputed values as data.

#' Write one uncertainty image per taxon
#'
#' @description
#' Every taxon gets a PNG under `uncertainty_png/`, at `dpi`.
#'
#' @param run The `mc_run` object.
#' @param out_dir Directory to write into.
#' @param dpi Resolution, in dots per inch.
#' @param say Function used to report progress.
#'
#' @return Character vector of the paths written, or `NULL` when there was
#'   nothing to draw.
#'
#' @keywords internal
#' @noRd
mc_write_uncertainty <- function(run, out_dir, dpi, say) {
    taxa <- unique(run$fit$pred_long$species)
    if (length(taxa) == 0) {
        return(NULL)
    }
    use_outliers <- isTRUE(run$fit$use_outliers)

    n_cells <- sum(run$fit$pred_long$species == taxa[1])
    mc_warn_page_count(length(taxa), n_cells)
    say("Drawing uncertainty for ", length(taxa), " taxa ...")

    # Working out a band refits the same FPCA models the run already fitted,
    # so the engine repeats the notes it made then. Those were reported in
    # aggregate at the end of the fit; repeating them per taxon here would
    # say nothing new and bury the console.
    named <- withCallingHandlers(
        {
            pages <- lapply(taxa, function(sp) {
                mc_taxon_pages(run$fit, sp, run$design, use_outliers)
            })
            stats::setNames(pages, taxa)
        },
        warning = function(w) invokeRestart("muffleWarning")
    )

    named <- mc_flatten_pages(named)
    if (length(named) == 0) {
        return(NULL)
    }

    mc_uncertainty_png(named, out_dir, dpi)
}

#' Flatten per-taxon page lists into one named list of pages
#'
#' @param pages Named list, one entry per taxon, each a list of plots.
#'
#' @return A named list of plots. A taxon with one imputed value keeps its
#'   own name; one with several gets a suffix per value.
#'
#' @keywords internal
#' @noRd
mc_flatten_pages <- function(pages) {
    out <- list()
    for (sp in names(pages)) {
        ps <- pages[[sp]]
        if (length(ps) == 0) next
        nms <- if (length(ps) == 1) {
            sp
        } else {
            paste0(sp, "_", seq_along(ps))
        }
        out[nms] <- ps
    }
    out
}

#' Write one PNG per taxon at the requested resolution
#'
#' @param named Named list of plots, one per taxon.
#' @param out_dir Directory to write into.
#' @param dpi Resolution in dots per inch.
#'
#' @return Character vector of the paths written.
#'
#' @keywords internal
#' @noRd
mc_uncertainty_png <- function(named, out_dir, dpi) {
    dir <- file.path(out_dir, "uncertainty_png")
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)

    files <- file.path(dir, paste0(mc_safe_name(names(named)), ".png"))
    for (i in seq_along(named)) {
        ggplot2::ggsave(
            files[i], named[[i]],
            width = 8, height = 5, dpi = dpi, units = "in"
        )
    }
    files
}

#' Make a taxon name usable as a file name
#'
#' @param x Character vector of taxon names.
#'
#' @return `x` with anything awkward for a file system replaced.
#'
#' @keywords internal
#' @noRd
mc_safe_name <- function(x) {
    out <- gsub("[^A-Za-z0-9._-]+", "_", x)
    out <- gsub("^_+|_+$", "", out)
    out[!nzchar(out)] <- "taxon"
    make.unique(out, sep = "_")
}

#' Check that a dpi argument is usable
#'
#' @param dpi As supplied by the caller.
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
#' @noRd
mc_check_dpi <- function(dpi) {
    ok <- is.numeric(dpi) && length(dpi) == 1 && !is.na(dpi) && dpi > 0
    if (!ok) {
        stop("dpi must be a single positive number.", call. = FALSE)
    }
    invisible(NULL)
}
