# Plotting one taxon.
#
# Everything drawn here is already in the run: mc_imputed_trajectory() gives
# one row per subject and time point, with the value, whether it was observed
# or imputed, and the interval for the imputed ones. So nothing is refitted
# and nothing is recomputed, which means the picture cannot disagree with the
# tables the run wrote.
#
# The earlier version of this function drew the fitted cluster of a single
# cell and spoke the fitting code's own vocabulary: internal subject codes
# and time positions. It is replaced because the question a reader actually
# has is simpler - where does this subject's trajectory go, which of its
# points were measured, and how sure is the model about the ones that were
# not - and because taking the caller's own subject names means the plot and
# the result tables can be talked about together.


#' Plot one taxon's trajectories, with the imputed values marked
#'
#' @description
#' Draws every subject's trajectory for one taxon, with one subject brought
#' to the front: its observed and imputed points are drawn with different
#' shapes, and the interval the run reported is shown on its imputed points
#' only.
#'
#' @details
#' The other subjects are drawn in grey, as context for the shape the model
#' had to work with. They carry no points and no intervals, because the
#' question the plot answers is about one subject.
#'
#' The interval is the one already in the run, the same numbers
#' [mc_uncertainty()] returns and `imputed_uncertainty.tsv` holds. It is a
#' prediction interval: it includes the residual variance, so it is wider
#' than a confidence interval for the mean curve would be.
#'
#' Intervals are shown for imputed points only, one bar each. An observed
#' value is a measurement, not an estimate, so the model has no uncertainty
#' to report about it, and a band joining the imputed points would be drawn
#' straight over the observed ones that lie between them.
#'
#' Nothing is fitted here, so no clustering and no outlier screening enters
#' into it. Whatever the run did, this draws what the run produced.
#'
#' @param run An `mc_run` object, as returned by [mc_run()].
#' @param taxon Character name of one taxon, as it appears in
#'   `run$completed`'s first column.
#' @param subject The subject to bring to the front, named as it is in
#'   your own metadata. Required: the plot is about one subject, and
#'   choosing one for the caller would make which subject they were looking
#'   at a matter of what happened to have the most gaps.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' demo <- mc_demo_data(n_taxa = 3, n_subjects = 6)
#' run <- suppressWarnings(mc_run(
#'     demo$abundance, demo$metadata,
#'     sample_col = "sample", subject_col = "subject", time_col = "time",
#'     C = 1, verbose = FALSE
#' ))
#'
#' tx <- rownames(demo$abundance)[1]
#'
#' # the subject named as your own metadata names it
#' subj <- run$imputed$subject[1]
#' mc_plot(run, tx, subject = subj)
#'
#' @seealso [mc_imputed_trajectory()] for the table behind the plot, and
#'   [mc_uncertainty()] for the intervals alone.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_errorbar
#'   labs theme_minimal scale_color_manual scale_shape_manual
#'   scale_x_continuous
#' @export
mc_plot <- function(run, taxon, subject) {
    if (missing(subject)) {
        stop(
            "subject is required: mc_plot() draws one subject's ",
            "trajectory. run$imputed$subject lists the subjects with ",
            "imputed values.",
            call. = FALSE
        )
    }
    if (!inherits(run, "mc_run")) {
        stop(
            "run must be an mc_run object, as returned by mc_run(). ",
            "Earlier versions of mc_plot() took a fit from mc_fit(); it ",
            "now takes the run, so that subjects and time points can be ",
            "named as your own metadata names them.",
            call. = FALSE
        )
    }

    tr <- mc_imputed_trajectory(run, taxon)
    subject <- mc_plot_subject(tr, subject)

    target <- tr[tr$subject == subject, , drop = FALSE]
    others <- tr[tr$subject != subject, , drop = FALSE]
    target <- target[order(target$time), , drop = FALSE]

    imp <- target[target$source == "imputed" & !is.na(target$lower), ,
        drop = FALSE
    ]

    p <- ggplot2::ggplot()

    # The other subjects first, so the one being looked at is drawn over
    # them rather than lost among them.
    if (nrow(others) > 0) {
        p <- p + ggplot2::geom_line(
            data = others[order(others$subject, others$time), ,
                drop = FALSE
            ],
            ggplot2::aes(
                x = .data$time, y = .data$value, group = .data$subject,
                colour = "other subjects"
            ),
            linewidth = 0.3, alpha = 0.45, na.rm = TRUE
        )
    }

    p <- p + ggplot2::geom_line(
        data = target,
        ggplot2::aes(
            x = .data$time, y = .data$value, colour = "this subject"
        ),
        linewidth = 0.9, na.rm = TRUE
    )

    # One bar per imputed point, not a band across them. A band would be
    # drawn over the observed points that lie between two imputed ones, and
    # so would claim uncertainty about values that were measured.
    if (nrow(imp) > 0) {
        p <- p + ggplot2::geom_errorbar(
            data = imp,
            ggplot2::aes(
                x = .data$time, ymin = .data$lower, ymax = .data$upper
            ),
            width = 0.12, linewidth = 0.5, colour = "#D55E00", na.rm = TRUE
        )
    }

    # Shape carries observed against imputed, so the two are told apart in
    # print and by a reader who cannot rely on colour.
    p <- p + ggplot2::geom_point(
        data = target,
        ggplot2::aes(
            x = .data$time, y = .data$value, shape = .data$source,
            colour = "this subject"
        ),
        size = 2.6, na.rm = TRUE
    )

    # The run works in positions along the time axis, but the caller may
    # have written the time points as labels. Where the two differ, the
    # labels are what the ticks say, so the picture reads in the caller's
    # own terms.
    ticks <- unique(tr[, c("time", "time_label")])
    ticks <- ticks[order(ticks$time), , drop = FALSE]
    if (!identical(as.character(ticks$time), as.character(ticks$time_label))) {
        p <- p + ggplot2::scale_x_continuous(
            breaks = ticks$time, labels = ticks$time_label
        )
    }

    n_sub <- length(unique(tr$subject))
    p +
        ggplot2::scale_color_manual(
            values = c(
                "this subject" = "#0072B2", "other subjects" = "grey65"
            ),
            breaks = c("this subject", "other subjects"),
            name = NULL
        ) +
        ggplot2::scale_shape_manual(
            values = c(observed = 16, imputed = 17),
            breaks = c("observed", "imputed"),
            name = NULL
        ) +
        ggplot2::labs(
            title = mc_plot_title(taxon),
            subtitle = paste0(
                "subject ", subject, " of ", n_sub, ", with ", nrow(imp),
                " imputed value(s) and their intervals"
            ),
            x = "time",
            y = "abundance"
        ) +
        ggplot2::theme_minimal()
}


#' Settle which subject the plot is about
#'
#' @param tr A trajectory table from [mc_imputed_trajectory()].
#' @param subject The caller's choice, or `NULL`.
#'
#' @return A single subject name.
#'
#' @keywords internal
#' @noRd
mc_plot_subject <- function(tr, subject) {
    subject <- as.character(subject)
    if (length(subject) != 1L) {
        stop("subject must be a single name.", call. = FALSE)
    }
    if (!subject %in% tr$subject) {
        stop(
            "There is no subject '", subject, "' in this run. It has ",
            mc_fmt_some(unique(tr$subject), n = 8), ".",
            call. = FALSE
        )
    }
    subject
}


#' Shorten a taxon name enough to sit in a title
#'
#' @description
#' A species-level lineage runs to well over a hundred characters, which no
#' title can hold. The last two ranks carry what a reader needs; the rest is
#' in the table.
#'
#' @param taxon The full taxon name.
#'
#' @return A single string.
#'
#' @keywords internal
#' @noRd
mc_plot_title <- function(taxon) {
    parts <- strsplit(taxon, "|", fixed = TRUE)[[1]]
    if (length(parts) <= 2L) {
        return(taxon)
    }
    paste(utils::tail(parts, 2), collapse = " | ")
}

