############################################
# Figures for mc_trends() and mc_compare().
#
# Drawing is a separate step, not something either function does: the point
# of a table of q-values is to be read before anything is drawn. mc_trends()
# therefore writes no figures at all, and mc_trend_plot() takes its result.
#
# The curves are population curves: the subject terms are excluded from the
# prediction, so what is drawn is the trajectory the model attributes to
# time and group, not to whichever subject happened to sit in the
# prediction grid.
#
# Models are not kept by either function, so each figure refits its own
# taxon. That is deterministic - mgcv::gam() with REML uses no random
# numbers - and far cheaper than storing the fits: a fitted gam at a few
# hundred samples is around half a megabyte, against seven kilobytes of
# coefficients that are actually wanted.
############################################


#' Refit one taxon, for a figure or a follow-up
#'
#' @param obj An `mc_trends` or `mc_compare` object.
#' @param taxon The taxon's name.
#'
#' @return A list of the fitted `gam` and the data it was fitted to.
#'
#' @keywords internal
#' @noRd
mc_trend_refit <- function(obj, taxon) {
    i <- match(taxon, rownames(obj$abundance))
    if (is.na(i)) {
        stop("No taxon '", taxon, "' in this result.", call. = FALSE)
    }
    s <- obj$settings

    d <- obj$design
    d$value <- as.numeric(obj$abundance[i, ])
    d <- d[!is.na(d$value), , drop = FALSE]
    d$subject <- factor(d$subject)
    if (!identical(s$group, "none")) d$mc_grp <- droplevels(d$mc_grp)

    # The fit behind the figure has to be the fit behind the p-value, so
    # the basis size comes from the row rather than from the settings:
    # with candidate k the settings hold the candidates, not the choice.
    kk <- obj$results$k[match(taxon, obj$results$taxon)]
    if (length(kk) != 1L || is.na(kk)) kk <- max(s$k)
    kk <- min(kk, length(unique(d$time)) - 1L)
    # The family too has to be the one behind the p-value, and with it the
    # offset: a count fit redrawn without its library sizes is a different
    # model from the one that was tested.
    fam_name <- s$family %||% "gaussian"
    use_off <- mc_trend_has_offset(fam_name)
    f <- mc_trend_formula(
        mc_trend_terms(
            kk, s$subject_effect, s$covariates, s$group, TRUE, use_off
        )
    )
    list(fit = mc_trend_gam(f, d, mc_trend_family_obj(fam_name)), data = d)
}


#' The smooth labels that belong to the subject, not to the trajectory
#'
#' @description
#' A population curve is wanted, so the subject terms are excluded from the
#' prediction. Their labels depend on `subject_effect`, so they are read off
#' the fit rather than assumed.
#'
#' @param m A fitted `gam`.
#'
#' @return A character vector of smooth labels, possibly empty.
#'
#' @keywords internal
#' @noRd
mc_trend_subject_smooths <- function(m) {
    if (length(m$smooth) == 0) {
        return(character(0))
    }
    labs <- vapply(m$smooth, function(s) s$label, character(1))
    labs[grepl("subject", labs, fixed = TRUE)]
}


#' A prediction grid over time, one row per group
#'
#' @param d The data the fit used.
#' @param group Whether there is a group term.
#' @param covariates Character vector, or `NULL`.
#' @param n_points How fine the grid is.
#'
#' @return A data.frame to pass to `predict()`.
#'
#' @keywords internal
#' @noRd
mc_trend_grid <- function(d, group, covariates = NULL, n_points = 200) {
    tt <- seq(min(d$time), max(d$time), length.out = n_points)
    nd <- if (identical(group, "none")) {
        data.frame(time = tt)
    } else {
        expand.grid(time = tt, mc_grp = levels(d$mc_grp))
    }
    # Any level will do: the subject terms are excluded from the prediction.
    nd$subject <- d$subject[1]
    # An offset cannot be excluded the way a smooth can, so the curve is
    # drawn at one stated library size rather than at none. The median is
    # the depth a typical sample was sequenced to, which makes the drawn
    # curve an expected count for such a sample.
    if (!is.null(d$mc_off)) {
        nd$mc_off <- stats::median(d$mc_off, na.rm = TRUE)
    }
    # A covariate has no shape of its own, so it is held at one value and
    # the drawn curve is the trajectory at that value.
    for (cv in covariates) {
        v <- d[[cv]]
        nd[[cv]] <- if (is.numeric(v)) stats::median(v, na.rm = TRUE) else {
            sort(unique(as.character(v)))[1]
        }
    }
    nd
}


#' One taxon's figure
#'
#' @description
#' The fitted trajectory with a 95 percent band and the observed points
#' behind it. A grouped comparison draws one curve per group.
#'
#' @param obj An `mc_trends` or `mc_compare` object.
#' @param taxon The taxon's name.
#'
#' @return A `ggplot`.
#'
#' @keywords internal
#' @noRd
mc_trend_panel <- function(obj, taxon) {
    s <- obj$settings
    rf <- mc_trend_refit(obj, taxon)
    m <- rf$fit
    d <- rf$data

    nd <- mc_trend_grid(d, s$group, s$covariates)
    pr <- stats::predict(
        m, nd,
        se.fit = TRUE, exclude = mc_trend_subject_smooths(m),
        newdata.guaranteed = FALSE
    )
    # The band is formed on the link scale, where the standard error is the
    # one mgcv reports and the interval is symmetric, and only then mapped
    # back. Under a log link that makes it asymmetric on the drawn scale,
    # which is correct: an interval symmetric around a count would run
    # negative. For the identity link of a gaussian fit the map is identity
    # and this is the arithmetic it has always been.
    li <- m$family$linkinv
    eta <- as.numeric(pr$fit)
    se <- as.numeric(pr$se.fit)
    nd$fit <- li(eta)
    nd$lo <- li(eta - 1.96 * se)
    nd$hi <- li(eta + 1.96 * se)

    if (identical(s$group, "none")) {
        # The colours are mapped to constants rather than set directly, so
        # that each mark earns a key. Without them a reader is left to
        # guess what the points are, and they are not all observations:
        # whatever table was fitted is what is drawn, so a completed table
        # puts its imputed cells here too.
        p <- ggplot2::ggplot(nd, ggplot2::aes(x = .data$time)) +
            ggplot2::geom_point(
                data = d,
                ggplot2::aes(y = .data$value, colour = "sample values"),
                alpha = 0.5, size = 1.3
            ) +
            ggplot2::geom_ribbon(
                ggplot2::aes(
                    ymin = .data$lo, ymax = .data$hi,
                    fill = mc_trend_band_label()
                ),
                alpha = 0.25
            ) +
            ggplot2::geom_line(
                ggplot2::aes(y = .data$fit, colour = "fitted trend"),
                linewidth = 1.1
            ) +
            ggplot2::scale_color_manual(
                values = c(
                    "sample values" = "black", "fitted trend" = "steelblue"
                ),
                breaks = c("fitted trend", "sample values"),
                name = NULL
            ) +
            ggplot2::scale_fill_manual(
                values = stats::setNames("steelblue",
                    mc_trend_band_label()),
                name = NULL
            ) +
            ggplot2::guides(
                colour = ggplot2::guide_legend(
                    override.aes = list(
                        linetype = c("solid", "blank"),
                        shape = c(NA, 16), alpha = 1, size = 2
                    )
                )
            ) +
            ggplot2::labs(
                title = mc_plot_title(taxon),
                x = "Time point", y = mc_trend_y_label(s)
            ) +
            mc_trend_theme()
        return(p)
    }

    curves <- ggplot2::ggplot(nd, ggplot2::aes(x = .data$time)) +
        ggplot2::geom_point(
            data = d,
            ggplot2::aes(y = .data$value, colour = .data$mc_grp),
            alpha = 0.35, size = 1.1
        ) +
        ggplot2::geom_ribbon(
            ggplot2::aes(
                ymin = .data$lo, ymax = .data$hi, fill = .data$mc_grp
            ),
            alpha = 0.22, colour = NA
        ) +
        ggplot2::geom_line(
            ggplot2::aes(y = .data$fit, colour = .data$mc_grp),
            linewidth = 0.9
        ) +
        ggplot2::labs(
            title = mc_plot_title(taxon),
            x = "Time point", y = mc_trend_y_label(s),
            colour = s$group_col, fill = s$group_col
        ) +
        mc_trend_theme()

    curves
}


#' Write one figure
#'
#' @param p The plot.
#' @param dir Directory to write into.
#' @param taxon The taxon's name, used for the file name.
#' @param dpi Resolution, in dots per inch.
#'
#' @return The path written.
#'
#' @keywords internal
#' @noRd
mc_trend_save <- function(p, dir, taxon, dpi) {
    f <- file.path(dir, paste0(mc_trend_stem(taxon), ".png"))
    grDevices::png(f, width = 7 * dpi, height = 6 * dpi, res = dpi)
    print(p)
    grDevices::dev.off()
    f
}


#' The band drawn around a fitted trend
#'
#' @description
#' A confidence band for the fitted curve, from `predict(se.fit = TRUE)`.
#' It is not the interval [mc_plot()] draws: that one covers a single
#' unobserved value and so carries the residual variance as well, which
#' makes it much wider. Naming them differently keeps the two pictures from
#' being read as the same thing.
#'
#' @return A single string.
#'
#' @keywords internal
#' @noRd
mc_trend_band_label <- function() "95% confidence band"


#' The theme the trend figures share
#'
#' @description
#' Centred title, and a base size that keeps a long lineage readable at the
#' sizes these are written out at.
#'
#' @return A ggplot2 theme.
#'
#' @keywords internal
#' @noRd
mc_trend_theme <- function() {
    ggplot2::theme_minimal(base_size = 13) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(
                hjust = 0.5, face = "bold", size = ggplot2::rel(1)
            ),
            plot.title.position = "plot"
        )
}


#' Draw the fitted trend of one taxon, or of several
#'
#' @description
#' Takes a finished [mc_trends()] or [mc_compare()] result and draws it.
#' Neither of those functions draws anything itself, so the table of
#' q-values can be read first and only the taxa worth seeing drawn.
#'
#' @details
#' The curve is a population curve: the subject terms are excluded from the
#' prediction, so what is drawn is the trajectory the model attributes to
#' time, and to group where there is one, rather than to whichever subject
#' sat in the prediction grid.
#'
#' The band is a 95% confidence band for that curve. It is not the interval
#' [mc_plot()] draws, which covers one unobserved value and is much wider.
#'
#' The points are every row that entered the fit. When the result was
#' produced from `mc_run()`'s completed table they include the imputed
#' cells, which is why they are labelled as sample values rather than as
#' observations.
#'
#' Neither function keeps its models, so this refits the one taxon. That is
#' deterministic, since `mgcv::gam()` with REML uses no random numbers.
#'
#' @param trend An `mc_trends` or `mc_compare` object.
#' @param taxa One taxon name, or a character vector of them, as they appear
#'   in `trend$results$taxon`. Each taxon gets a figure of its own; nothing
#'   is combined into a grid, since a panel titled with a lineage is
#'   unreadable at a quarter of the page. `NULL`, the default, takes the
#'   taxa that reached `q < alpha`. A name is drawn whatever its q-value:
#'   asking for one by name is already a decision that it is worth seeing.
#' @param out_dir Directory to write the figures into, one PNG per taxon,
#'   or `NULL` to return them without writing. Created if it does not exist.
#' @param dpi Resolution of the written images, in dots per inch.
#' @param verbose Logical. Whether to report what was written.
#'
#' @return One taxon gives the figure itself, a `ggplot`. Several give a
#'   list of them, named by taxon, so `p[[2]]` or `print(p)` draws them one
#'   at a time. When `out_dir` was written to, the return is invisible.
#'
#' @examples
#' demo <- mc_demo_data(n_taxa = 3, n_subjects = 8)
#' tr <- mc_trends(
#'     demo$abundance, demo$metadata,
#'     sample_col = "sample", subject_col = "subject", time_col = "time",
#'     verbose = FALSE
#' )
#'
#' # one taxon, drawn on the spot
#' mc_trend_plot(tr, tr$results$taxon[1])
#'
#' # several, one figure each
#' ps <- mc_trend_plot(tr, tr$results$taxon[1:2])
#' ps[[1]]
#'
#' @seealso [mc_trends()] and [mc_compare()] to produce the result,
#'   [mc_plot()] for one subject's own trajectory.
#'
#' @export
mc_trend_plot <- function(trend, taxa = NULL, out_dir = NULL, dpi = 300,
                          verbose = TRUE) {
    if (!inherits(trend, c("mc_trends", "mc_compare"))) {
        stop(
            "trend must be an mc_trends or mc_compare object, as returned ",
            "by mc_trends() or mc_compare(). For one subject's trajectory ",
            "from a run, use mc_plot().",
            call. = FALSE
        )
    }
    want <- mc_trend_which(trend, taxa)

    if (!is.null(out_dir)) {
        mc_check_dpi(dpi)
        if (!dir.exists(out_dir)) {
            dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        }
    }

    ps <- stats::setNames(vector("list", length(want)), want)
    files <- character(0)
    for (i in seq_along(want)) {
        # One taxon lets its own error through, since that error is the
        # answer. A batch does not: losing forty-nine figures to the
        # fiftieth would be the worse failure, so it warns and carries on.
        ps[[i]] <- if (length(want) == 1L) {
            mc_trend_panel(trend, want[i])
        } else {
            tryCatch(mc_trend_panel(trend, want[i]), error = function(e) {
                warning(
                    "could not draw '", want[i], "': ", conditionMessage(e),
                    call. = FALSE
                )
                NULL
            })
        }
        if (!is.null(out_dir) && !is.null(ps[[i]])) {
            files <- c(files, mc_trend_save(ps[[i]], out_dir, want[i], dpi))
        }
    }
    ps <- ps[!vapply(ps, is.null, logical(1))]

    if (isTRUE(verbose) && length(files) > 0) {
        message("  wrote ", length(files), " figure(s) to ", out_dir)
    }

    out <- if (length(ps) == 1L) ps[[1]] else ps
    if (is.null(out_dir)) out else invisible(out)
}


#' Which taxa a figure was asked for
#'
#' @param trend The result object.
#' @param taxa What the caller passed, possibly `NULL`.
#'
#' @return A character vector of taxon names, each one drawable.
#'
#' @keywords internal
#' @noRd
mc_trend_which <- function(trend, taxa) {
    r <- trend$results
    if (is.null(taxa)) {
        hit <- r$taxon[!is.na(r$q_value) & r$q_value < trend$settings$alpha]
        if (length(hit) == 0) {
            stop(
                "No taxon reached q < ", trend$settings$alpha,
                ", so there is nothing to draw by default. Name the taxa ",
                "to draw, for instance ",
                "mc_trend_plot(trend, trend$results$taxon[1]).",
                call. = FALSE
            )
        }
        return(hit)
    }

    if (!is.character(taxa) || length(taxa) == 0) {
        stop(
            "taxa must be one taxon name or a character vector of them.",
            call. = FALSE
        )
    }
    taxa <- unique(taxa)

    gone <- setdiff(taxa, r$taxon)
    if (length(gone) > 0) {
        stop(
            "Not in this result: ", mc_fmt_some(gone), ". It has ",
            nrow(r), " taxa, among them ",
            mc_fmt_some(r$taxon, n = 4), ".",
            call. = FALSE
        )
    }

    # A taxon that never fitted has no fit to draw, and refitting it here
    # would fail for the reason already recorded in its note.
    ok <- r$converged[match(taxa, r$taxon)]
    bad <- taxa[is.na(ok) | !ok]
    if (length(bad) > 0) {
        note <- r$note[match(bad[1], r$taxon)]
        stop(
            if (length(bad) == 1L) {
                paste0("Taxon '", bad, "' was not tested")
            } else {
                paste0("Not tested: ", mc_fmt_some(bad))
            },
            if (!is.na(note)) paste0(" (", note, ")") else "",
            ".",
            call. = FALSE
        )
    }
    taxa
}


#' A file name for one taxon
#'
#' @description
#' A full lineage runs well past any sensible file name, and what
#' distinguishes two of them is at the end: cutting the first 80 characters
#' gives two species of one genus the same name, and the second figure then
#' overwrites the first. So the name is built from the end of the lineage,
#' and a short hash of the whole string keeps two taxa apart even when
#' their last ranks read alike.
#'
#' @param taxon The full taxon name.
#'
#' @return A single string, safe as a file name on every platform.
#'
#' @keywords internal
#' @noRd
mc_trend_stem <- function(taxon) {
    parts <- strsplit(taxon, "|", fixed = TRUE)[[1]]
    tail_part <- paste(utils::tail(parts, 2), collapse = "_")

    stem <- gsub("[^A-Za-z0-9._-]+", "_", tail_part)
    stem <- gsub("^_+|_+$", "", stem)
    stem <- substr(stem, 1, 60)
    if (!nzchar(stem)) stem <- "taxon"

    # Eight hex digits of the whole lineage. Short enough to read past,
    # long enough that a collision is not something to plan for.
    paste0(stem, "_", substr(mc_trend_hash(taxon), 1, 8))
}


#' A stable short hash of a string
#'
#' @description
#' Used only to keep file names apart, so it has to be deterministic
#' between sessions and nothing more. A polynomial rolling hash, which
#' needs no package.
#'
#' @details
#' Deliberately arithmetic rather than bitwise. R has no unsigned integer,
#' and `bitwXor()` returns `NA` for a value past 2^31, so an FNV-style hash
#' written the obvious way silently produces NA file names. Here every
#' intermediate stays below 2^31 * 31 + 255, which a double holds exactly.
#'
#' @param x A single string.
#'
#' @return An eight-character hex string.
#'
#' @keywords internal
#' @noRd
mc_trend_hash <- function(x) {
    h <- 0
    for (b in as.integer(charToRaw(x))) {
        h <- (h * 31 + b) %% 2147483647
    }
    sprintf("%08x", as.integer(h))
}
