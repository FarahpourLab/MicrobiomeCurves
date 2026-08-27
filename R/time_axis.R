# Turning a metadata time column into a numeric axis.
#
# The model needs a number for each time point, because it fits a curve over
# time. Studies do not always record one: a column may hold "baseline",
# "week1", "week4", or a factor with its own level order.
#
# Numbers are used as they stand, so real spacing is respected. Labels are
# placed at consecutive positions in the order they are given, which assumes
# equal spacing between them. That assumption changes the shape of every
# fitted curve, so it is stated out loud rather than made quietly.

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
tti_time_axis <- function(time_raw, time_col) {
    if (anyNA(time_raw)) {
        stop(
            "metadata column '", time_col, "' contains missing time points. ",
            "Every sample needs a time.",
            call. = FALSE
        )
    }

    if (is.numeric(time_raw)) {
        return(tti_axis_numeric(as.numeric(time_raw)))
    }

    txt <- trimws(as.character(time_raw))
    if (all(tti_looks_numeric(txt)) && !is.factor(time_raw)) {
        return(tti_axis_numeric(as.numeric(txt)))
    }

    tti_axis_labels(time_raw, txt, time_col)
}

#' Use a genuinely numeric time column as it stands
#'
#' @param time Numeric vector of time points.
#'
#' @return A time axis list, with `literal = TRUE`.
#'
#' @keywords internal
#' @noRd
tti_axis_numeric <- function(time) {
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
tti_axis_labels <- function(time_raw, txt, time_col) {
    lv <- if (is.factor(time_raw)) {
        levels(time_raw)[levels(time_raw) %in% txt]
    } else {
        unique(txt)
    }

    if (length(lv) < 2) {
        stop(
            "metadata column '", time_col, "' has only one distinct value (",
            tti_fmt_some(lv), "). Imputation over time needs at least two.",
            call. = FALSE
        )
    }

    shown <- tti_fmt_some(paste0(lv, " = ", seq_along(lv)), n = 8)
    source_of_order <- if (is.factor(time_raw)) {
        "the factor's levels"
    } else {
        "the order the rows appear in"
    }

    warning(
        "metadata column '", time_col, "' is not numeric, so its values are ",
        "placed in order at equal spacing: ", shown,
        ". The order comes from ", source_of_order,
        ". If the real intervals are uneven, supply them as numbers instead, ",
        "since the fitted curves depend on the spacing.",
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
#' @param axis A time axis list from [tti_time_axis()].
#' @param time Numeric vector of positions on that axis.
#'
#' @return Character vector of labels, one per entry of `time`.
#'
#' @keywords internal
#' @noRd
tti_time_label <- function(axis, time) {
    if (is.null(axis) || isTRUE(axis$literal)) {
        return(as.character(time))
    }
    idx <- match(time, axis$positions)
    ifelse(is.na(idx), as.character(time), axis$levels[idx])
}

#' Describe a time axis in one line
#'
#' @param axis A time axis list from [tti_time_axis()].
#'
#' @return A single string naming the time points in order.
#'
#' @keywords internal
#' @noRd
tti_axis_text <- function(axis) {
    if (isTRUE(axis$literal)) {
        return(tti_fmt_some(axis$positions, n = 12))
    }
    tti_fmt_some(
        paste0(axis$levels, " (", axis$positions, ")"),
        n = 12
    )
}
