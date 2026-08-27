# One page per taxon, its imputed values side by side.
#
# The band differs from value to value, because each is the prediction
# interval for a different subject, so the values are drawn as facets rather
# than overlaid. A taxon with more imputed values than fit legibly on one
# page spills onto further pages rather than having any of them dropped.

# Facets per page. Six fits a landscape page at a readable size; beyond that
# the panels are too small to read a band off.
TTI_FACETS_PER_PAGE <- 6L

#' Every page for one taxon
#'
#' @param fit The `tti_fit`.
#' @param species_name Character name of the taxon.
#' @param design The `tti_design`.
#' @param use_outliers Logical. Whether the fit screened for outliers.
#'
#' @return A list of `ggplot` objects: one page when the taxon's imputed
#'   values fit, more when they do not.
#'
#' @keywords internal
#' @noRd
tti_taxon_pages <- function(fit, species_name, design, use_outliers) {
    pl <- fit$pred_long
    cells <- pl[pl$species == species_name, ]
    if (nrow(cells) == 0) {
        return(list())
    }

    panels <- lapply(seq_len(nrow(cells)), function(i) {
        tti_panel_data(
            fit, species_name, cells$rep[i], cells$time[i],
            design, use_outliers
        )
    })
    panels <- Filter(Negate(is.null), panels)
    if (length(panels) == 0) {
        return(list())
    }

    groups <- split(
        panels,
        ceiling(seq_along(panels) / TTI_FACETS_PER_PAGE)
    )
    n_pages <- length(groups)

    lapply(seq_along(groups), function(i) {
        tti_faceted_page(
            groups[[i]], species_name, use_outliers,
            if (n_pages > 1) paste0(" (page ", i, " of ", n_pages, ")") else ""
        )
    })
}

#' Draw one page of faceted panels
#'
#' @param panels List of panel data, one per imputed value.
#' @param species_name Character name of the taxon.
#' @param use_outliers Logical. Whether the fit screened for outliers.
#' @param suffix Character added to the title when a taxon spans pages.
#'
#' @return A `ggplot` object.
#'
#' @keywords internal
#' @noRd
tti_faceted_page <- function(panels, species_name, use_outliers, suffix) {
    merged <- tti_merge_panels(panels)

    p <- tti_uncertainty_panel(merged, use_outliers) +
        ggplot2::facet_wrap(~cell, scales = "free_y") +
        ggplot2::labs(
            title = paste0(species_name, suffix),
            subtitle = panels[[1]]$subtitle
        )

    if (length(panels) > 1) {
        # One dashed marker per facet, each at its own target time.
        p <- p + ggplot2::geom_vline(
            data = merged$targets,
            ggplot2::aes(xintercept = .data$target_time),
            linewidth = 0.3, linetype = "22", colour = "grey40"
        )
    }
    p
}

#' Stack several panels' data, tagged by facet
#'
#' @param panels List of panel data from [tti_panel_data()].
#'
#' @return A single panel-data list, with a `cell` column throughout and a
#'   `targets` table giving each facet's target time.
#'
#' @keywords internal
#' @noRd
tti_merge_panels <- function(panels) {
    label <- vapply(panels, function(x) x$cell_label, character(1))
    label <- factor(label, levels = label)

    tag <- function(field) {
        parts <- lapply(seq_along(panels), function(i) {
            d <- panels[[i]][[field]]
            if (is.null(d) || nrow(d) == 0) {
                return(NULL)
            }
            d$cell <- label[i]
            d
        })
        out <- do.call(rbind, Filter(Negate(is.null), parts))
        if (is.null(out)) {
            data.frame(cell = factor(character(0), levels = levels(label)))
        } else {
            out
        }
    }

    first <- panels[[1]]
    list(
        bands = tag("bands"), traj = tag("traj"),
        observed = tag("observed"), imputed = tag("imputed"),
        truth = tag("truth"),
        targets = data.frame(
            cell = label,
            target_time = vapply(panels, function(x) x$target_time, numeric(1))
        ),
        # A single-facet page keeps the plain vline; multi-facet pages draw
        # their own, one per target.
        target_time = if (length(panels) == 1) {
            first$target_time
        } else {
            numeric(0)
        },
        title = first$title, subtitle = first$subtitle,
        colour_breaks = first$colour_breaks, fill_breaks = first$fill_breaks
    )
}

#' Warn when a run is about to draw a great many pages
#'
#' @param n_taxa Integer number of taxa.
#' @param n_cells Integer number of imputed values per taxon.
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
#' @noRd
tti_warn_page_count <- function(n_taxa, n_cells) {
    per_taxon <- ceiling(n_cells / TTI_FACETS_PER_PAGE)
    total <- n_taxa * per_taxon
    if (total <= 200) {
        return(invisible(NULL))
    }

    warning(
        "Drawing ", total, " uncertainty pages (", n_taxa, " taxa x ",
        per_taxon, " page(s) each, for ", n_cells,
        " imputed value(s) per taxon). This dominates the run time. Pass ",
        "plots = FALSE to skip it, or use tti_uncertainty() to read the ",
        "intervals for the taxa you care about.",
        call. = FALSE
    )
    invisible(NULL)
}
