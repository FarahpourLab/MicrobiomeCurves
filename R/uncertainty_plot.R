# Drawing per-taxon uncertainty, and collecting it into one PDF.
#
# One page per taxon, each imputed value shown with its interval. A page is
# readable on its own, and the file as a whole is the thing to look at before
# treating imputed values as data.

#' Draw the uncertainty of one taxon's imputed values
#'
#' @param unc Data frame from [tti_taxon_uncertainty()].
#' @param species_name Character name of the taxon, used in the title.
#'
#' @return A `ggplot` object.
#'
#' @keywords internal
#' @noRd
tti_plot_uncertainty <- function(unc, species_name) {
    unc$label <- paste0(unc$subject, " @ ", unc$time_label)
    unc <- unc[order(unc$subject, unc$time, method = "radix"), ]
    unc$label <- factor(unc$label, levels = unc$label)

    n_wide <- sum(!is.na(unc$se))
    sub <- if (n_wide == 0) {
        "no interval could be computed for these cells"
    } else {
        paste0(
            nrow(unc), " imputed value(s); bars are 95% analytic intervals"
        )
    }

    ggplot2::ggplot(
        unc,
        ggplot2::aes(x = .data$label, y = .data$imputed)
    ) +
        ggplot2::geom_errorbar(
            ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
            width = 0.2, na.rm = TRUE
        ) +
        ggplot2::geom_point(size = 2, na.rm = TRUE) +
        ggplot2::labs(
            title = species_name, subtitle = sub,
            x = NULL, y = "imputed value (CLR)"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )
}

#' Write one uncertainty page per taxon
#'
#' @description
#' Every taxon gets a page in one PDF. The PDF is vector art, so it is sharp
#' at any magnification and `dpi` does not apply to it; `dpi` governs the
#' PNG copies, which are written one per taxon when asked for.
#'
#' @param run The `tti_run` object.
#' @param out_dir Directory to write into.
#' @param dpi Resolution for the PNG copies, in dots per inch.
#' @param format Either `"pdf"`, `"png"`, or `"both"`.
#' @param say Function used to report progress.
#'
#' @return Character vector of the paths written, or `NULL` when there was
#'   nothing to draw.
#'
#' @keywords internal
#' @noRd
tti_write_uncertainty <- function(run, out_dir, dpi, format, say) {
    taxa <- unique(run$fit$pred_long$species)
    if (length(taxa) == 0) {
        return(NULL)
    }

    say("Drawing uncertainty for ", length(taxa), " taxa ...")

    # Working out an interval refits the same FPCA models the run already
    # fitted, so the engine repeats the notes it made then. Those were
    # reported in aggregate at the end of the fit; repeating them once per
    # taxon here would say nothing new and bury the console.
    plots <- withCallingHandlers(
        lapply(taxa, function(sp) {
            unc <- tti_taxon_uncertainty(run$fit, sp, run$design)
            if (is.null(unc)) NULL else tti_plot_uncertainty(unc, sp)
        }),
        warning = function(w) invokeRestart("muffleWarning")
    )
    named <- stats::setNames(plots, taxa)
    named <- named[!vapply(named, is.null, logical(1))]
    if (length(named) == 0) {
        return(NULL)
    }

    out <- character(0)
    if (format %in% c("pdf", "both")) {
        out <- c(out, tti_uncertainty_pdf(named, out_dir))
    }
    if (format %in% c("png", "both")) {
        out <- c(out, tti_uncertainty_png(named, out_dir, dpi))
    }
    out
}

#' Collect the pages into one PDF
#'
#' @param named Named list of plots, one per taxon.
#' @param out_dir Directory to write into.
#'
#' @return The path written.
#'
#' @keywords internal
#' @noRd
tti_uncertainty_pdf <- function(named, out_dir) {
    path <- file.path(out_dir, "uncertainty_by_taxon.pdf")

    grDevices::pdf(path, width = 8, height = 5, onefile = TRUE)
    on.exit(grDevices::dev.off(), add = TRUE)

    # Drawn through grid rather than print(), which draws the same page but
    # is the printing method for interactive use.
    for (p in named) {
        grid::grid.newpage()
        grid::grid.draw(ggplot2::ggplotGrob(p))
    }

    path
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
tti_uncertainty_png <- function(named, out_dir, dpi) {
    dir <- file.path(out_dir, "uncertainty_png")
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)

    files <- file.path(dir, paste0(tti_safe_name(names(named)), ".png"))
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
tti_safe_name <- function(x) {
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
tti_check_dpi <- function(dpi) {
    ok <- is.numeric(dpi) && length(dpi) == 1 && !is.na(dpi) && dpi > 0
    if (!ok) {
        stop("dpi must be a single positive number.", call. = FALSE)
    }
    invisible(NULL)
}
