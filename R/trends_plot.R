############################################
# Figures for mc_trends() and mc_compare().
#
# One figure per taxon that reached significance. The curves are population
# curves: the subject terms are excluded from the prediction, so what is
# drawn is the trajectory the model attributes to time and group, not to
# whichever subject happened to sit in the prediction grid.
#
# Models are not kept by either function, so each figure refits its own
# taxon. That is deterministic - mgcv::gam() with REML uses no random
# numbers - and far cheaper than storing the fits: a fitted gam at a few
# hundred samples is around half a megabyte, against seven kilobytes of
# coefficients that are actually wanted.
############################################


#' Draw a figure per significant taxon
#'
#' @param obj An `mc_trends` or `mc_compare` object.
#' @param out_dir Directory to write into. Figures go in a `plots`
#'   subdirectory of it.
#' @param dpi Resolution for PNG output.
#' @param plot_format One of `"pdf"`, `"png"` or `"both"`.
#' @param verbose Logical.
#'
#' @return Character vector of the paths written.
#'
#' @keywords internal
#' @noRd
mc_trend_plots <- function(obj, out_dir, dpi = 300, plot_format = "pdf",
                           verbose = TRUE) {
    r <- obj$results
    hit <- r$taxon[!is.na(r$q_value) & r$q_value < obj$settings$alpha]
    if (length(hit) == 0) {
        if (isTRUE(verbose)) {
            message("  no taxon reached q < ", obj$settings$alpha,
                    "; no figures drawn")
        }
        return(character(0))
    }

    dir <- file.path(out_dir, "plots")
    if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }

    paths <- character(0)
    for (tx in hit) {
        p <- tryCatch(mc_trend_panel(obj, tx), error = function(e) NULL)
        if (is.null(p)) next
        paths <- c(paths, mc_trend_save(p, dir, tx, dpi, plot_format))
    }
    if (isTRUE(verbose)) {
        message("  drew ", length(hit), " figure(s)")
    }
    paths
}


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

    kk <- min(s$k, length(unique(d$time)) - 1L)
    f <- mc_trend_formula(
        mc_trend_terms(kk, s$subject_effect, s$covariates, s$group, TRUE)
    )
    list(fit = mc_trend_gam(f, d), data = d)
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
#' The fitted trajectory with a 95 percent band, the observed points behind
#' it, and for a two-group comparison a second panel holding the difference
#' between the curves with the band that says where it leaves zero.
#'
#' @param obj An `mc_trends` or `mc_compare` object.
#' @param taxon The taxon's name.
#'
#' @return A `ggplot` or a `patchwork` of two.
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
    nd$fit <- as.numeric(pr$fit)
    nd$lo <- nd$fit - 1.96 * as.numeric(pr$se.fit)
    nd$hi <- nd$fit + 1.96 * as.numeric(pr$se.fit)

    row <- obj$results[obj$results$taxon == taxon, , drop = FALSE][1, ]
    sub <- sprintf(
        "edf %.2f   p %.3g   q %.3g   %d obs, %d subjects",
        row$edf, row$p_value, row$q_value, row$n_obs, row$n_subjects
    )

    if (identical(s$group, "none")) {
        p <- ggplot2::ggplot(nd, ggplot2::aes(x = .data$time)) +
            ggplot2::geom_point(
                data = d, ggplot2::aes(y = .data$value),
                colour = "grey45", alpha = 0.45, size = 1.1
            ) +
            ggplot2::geom_ribbon(
                ggplot2::aes(ymin = .data$lo, ymax = .data$hi),
                fill = "steelblue", alpha = 0.25
            ) +
            ggplot2::geom_line(
                ggplot2::aes(y = .data$fit),
                colour = "steelblue", linewidth = 0.9
            ) +
            ggplot2::labs(
                title = taxon, subtitle = sub,
                x = s$time_col, y = "abundance (CLR)"
            ) +
            ggplot2::theme_minimal()
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
            title = taxon, subtitle = sub,
            x = s$time_col, y = "abundance (CLR)",
            colour = s$group_col, fill = s$group_col
        ) +
        ggplot2::theme_minimal()

    diff <- mc_trend_diff_panel(obj, m, d)
    if (is.null(diff)) {
        return(curves)
    }
    patchwork::wrap_plots(curves, diff, ncol = 1, heights = c(2, 1))
}


#' The difference panel, for a two-group comparison
#'
#' @description
#' Only the ordered parameterisation holds a single difference smooth, so
#' only then is there one curve to draw. With `reference = "global"` there
#' are as many departures as there are groups and the comparison between any
#' two of them belongs in a pairwise follow-up, not here.
#'
#' @param obj The result object.
#' @param m The refitted model.
#' @param d The data it was fitted to.
#'
#' @return A `ggplot`, or `NULL` when there is no single difference term.
#'
#' @keywords internal
#' @noRd
mc_trend_diff_panel <- function(obj, m, d) {
    s <- obj$settings
    if (!identical(s$group, "ordered") || nlevels(d$mc_grp) != 2L) {
        return(NULL)
    }
    labs <- vapply(m$smooth, function(x) x$label, character(1))
    want <- labs[grepl("):", labs, fixed = TRUE)]
    if (length(want) != 1L) {
        return(NULL)
    }

    nd <- mc_trend_grid(d, s$group, s$covariates)
    nd <- nd[nd$mc_grp == levels(d$mc_grp)[2], , drop = FALSE]
    tr <- stats::predict(m, nd, type = "terms", terms = want, se.fit = TRUE)

    g <- data.frame(
        time = nd$time,
        diff = as.numeric(tr$fit[, 1]),
        se = as.numeric(tr$se.fit[, 1])
    )
    g$lo <- g$diff - 1.96 * g$se
    g$hi <- g$diff + 1.96 * g$se
    # Where the band clears zero is where the two trajectories are
    # distinguishable, which is the "when" that a single p-value cannot say.
    g$out <- g$lo > 0 | g$hi < 0

    ggplot2::ggplot(g, ggplot2::aes(x = .data$time)) +
        ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
        ggplot2::geom_ribbon(
            ggplot2::aes(ymin = .data$lo, ymax = .data$hi),
            fill = "grey30", alpha = 0.2
        ) +
        ggplot2::geom_line(
            ggplot2::aes(y = .data$diff),
            colour = "grey15", linewidth = 0.8
        ) +
        ggplot2::geom_point(
            data = g[g$out, , drop = FALSE],
            ggplot2::aes(y = .data$diff), colour = "firebrick", size = 0.7
        ) +
        ggplot2::labs(
            x = s$time_col,
            y = paste0(levels(d$mc_grp)[2], " - ", levels(d$mc_grp)[1]),
            caption = paste0(
                "red where the 95% band excludes zero; ",
                "covariates held at a reference value"
            )
        ) +
        ggplot2::theme_minimal()
}


#' Write one figure in the requested format(s)
#'
#' @param p The plot.
#' @param dir Directory to write into.
#' @param taxon The taxon's name, used for the file name.
#' @param dpi Resolution for PNG output.
#' @param plot_format One of `"pdf"`, `"png"` or `"both"`.
#'
#' @return Character vector of the paths written.
#'
#' @keywords internal
#' @noRd
mc_trend_save <- function(p, dir, taxon, dpi, plot_format) {
    # A taxonomy string carries slashes, semicolons and spaces, none of
    # which belong in a file name.
    stem <- gsub("[^A-Za-z0-9._-]+", "_", taxon)
    stem <- substr(sub("^_+|_+$", "", stem), 1, 80)
    if (!nzchar(stem)) stem <- "taxon"

    paths <- character(0)
    if (plot_format %in% c("pdf", "both")) {
        f <- file.path(dir, paste0(stem, ".pdf"))
        grDevices::pdf(f, width = 7, height = 6)
        print(p)
        grDevices::dev.off()
        paths <- c(paths, f)
    }
    if (plot_format %in% c("png", "both")) {
        f <- file.path(dir, paste0(stem, ".png"))
        grDevices::png(f, width = 7 * dpi, height = 6 * dpi, res = dpi)
        print(p)
        grDevices::dev.off()
        paths <- c(paths, f)
    }
    paths
}
