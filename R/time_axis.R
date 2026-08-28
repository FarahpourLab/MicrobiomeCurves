# Turning a metadata time column into a numeric axis.
#
# The model needs a number for each time point, because it fits a curve over
# time. Studies do not always record one: a column may hold "baseline",
# "week1", "week4", or a factor with its own level order.
#
# What the model actually receives is the ORDER of the time points, not
# their values: mc_encode_samples() replaces each time by its rank, so
# times 0, 1, 2, 60 are fitted exactly as 0, 1, 2, 3 would be. The values
# are kept here so reports, axes and column names can speak in the study's
# own terms, but they do not enter the fit.
#
# Labels are therefore no worse off than numbers: both are reduced to an
# order. What a label loses is only the ability to say what the intervals
# were, which the fit was not using anyway.

#' Build a numeric time axis from a metadata column
#'
#' @param time_raw The time column as stored: numeric, character, or factor.
#' @param time_col Its name, used in messages.
#'
#' @return List with `time`, the numeric position of each sample; `levels`,
#'   the ordered labels; `positions`, their numeric values; and `literal`,
#'   `TRUE` when the column held real numbers.
#'
#' @keywords internal
#' @noRd
mc_time_axis <- function(time_raw, time_col) {
    if (anyNA(time_raw)) {
        stop(
            "metadata column '", time_col, "' contains missing time points. ",
            "Every sample needs a time.",
            call. = FALSE
        )
    }

    if (is.numeric(time_raw)) {
        return(mc_axis_numeric(as.numeric(time_raw)))
    }

    txt <- trimws(as.character(time_raw))
    if (all(mc_looks_numeric(txt)) && !is.factor(time_raw)) {
        return(mc_axis_numeric(as.numeric(txt)))
    }

    mc_axis_labels(time_raw, txt, time_col)
}

#' Use a genuinely numeric time column as it stands
#'
#' @param time Numeric vector of time points.
#'
#' @return A time axis list, with `literal = TRUE`.
#'
#' @keywords internal
#' @noRd
mc_axis_numeric <- function(time) {
    lv <- sort(unique(time))
    list(
        time = time,
        levels = as.character(lv),
        positions = lv,
        literal = TRUE
    )
}

#' Place ordered labels on consecutive positions
#'
#' @description
#' A factor's level order is taken as the intended order. Anything else is
#' ordered by first appearance, which is the order the rows were written in.
#'
#' @param time_raw The column as stored, used for factor levels.
#' @param txt The column as character.
#' @param time_col Its name, used in messages.
#'
#' @return A time axis list, with `literal = FALSE`.
#'
#' @keywords internal
#' @noRd
mc_axis_labels <- function(time_raw, txt, time_col) {
    lv <- if (is.factor(time_raw)) {
        levels(time_raw)[levels(time_raw) %in% txt]
    } else {
        unique(txt)
    }

    if (length(lv) < 2) {
        stop(
            "metadata column '", time_col, "' has only one distinct value (",
            mc_fmt_some(lv), "). Imputation over time needs at least two.",
            call. = FALSE
        )
    }

    shown <- mc_fmt_some(paste0(lv, " = ", seq_along(lv)), n = 8)
    source_of_order <- if (is.factor(time_raw)) {
        "the factor's levels"
    } else {
        "the order the rows appear in"
    }

    warning(
        "metadata column '", time_col, "' is not numeric, so its values are ",
        "read as an order: ", shown,
        ". The order comes from ", source_of_order,
        ". Check it is the order you intend. Note that the fit uses the ",
        "order of time points rather than their values, so numeric times ",
        "would be treated the same way.",
        call. = FALSE
    )

    list(
        time = match(txt, lv),
        levels = lv,
        positions = seq_along(lv),
        literal = FALSE
    )
}

#' Name a time point the way the user wrote it
#'
#' @description
#' Positions are what the model works in, but a study that supplied labels
#' should see its labels back, in every report, column name and log line.
#'
#' @param axis A time axis list from [mc_time_axis()].
#' @param time Numeric vector of positions on that axis.
#'
#' @return Character vector of labels, one per entry of `time`.
#'
#' @keywords internal
#' @noRd
mc_time_label <- function(axis, time) {
    if (is.null(axis) || isTRUE(axis$literal)) {
        return(as.character(time))
    }
    idx <- match(time, axis$positions)
    ifelse(is.na(idx), as.character(time), axis$levels[idx])
}

#' Describe a time axis in one line
#'
#' @param axis A time axis list from [mc_time_axis()].
#'
#' @return A single string naming the time points in order.
#'
#' @keywords internal
#' @noRd
mc_axis_text <- function(axis) {
    if (isTRUE(axis$literal)) {
        return(mc_fmt_some(axis$positions, n = 12))
    }
    mc_fmt_some(
        paste0(axis$levels, " (", axis$positions, ")"),
        n = 12
    )
}
