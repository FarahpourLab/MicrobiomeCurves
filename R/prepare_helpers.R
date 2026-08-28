#' Default parser for `<subject>.<time>` column names
#'
#' @param cols Character vector of column names.
#'
#' @return Tibble with `col`, `rep` and `time` columns. `rep` and `time` are
#'   `NA` for names that do not match.
#'
#' @keywords internal
#' @noRd
mc_parse_cols <- function(cols) {
    parts <- stringr::str_match(cols, "^([^\\.]+)\\.(\\d+)$")
    tibble(
        col = cols,
        rep = parts[, 2],
        time = as.numeric(parts[, 3])
    )
}

#' Turn a mask list into subject-timepoint pairs
#'
#' @param mask_list Data frame whose first two columns are subject and time,
#'   or a list of length-2 vectors.
#'
#' @return Data frame with `rep` and `time` columns.
#'
#' @keywords internal
#' @noRd
mc_mask_pairs_from_list <- function(mask_list) {
    if (is.data.frame(mask_list)) {
        if (ncol(mask_list) < 2) {
            stop(
                "mask_list has ", ncol(mask_list), " column(s). As a data ",
                "frame it needs at least two, read as subject and time, ",
                "for example data.frame(rep = \"A1\", time = 3).",
                call. = FALSE
            )
        }
        ml <- mask_list[, seq_len(2)]
        names(ml) <- c("rep", "time")
        return(ml)
    }

    if (!is.list(mask_list)) {
        stop(
            "mask_list must be a data frame or a list, not ",
            class(mask_list)[1], ".",
            call. = FALSE
        )
    }

    mc_check_pair_shape(mask_list)

    tibble(
        rep = vapply(
            mask_list, function(x) as.character(x[1]), character(1)
        ),
        time = mc_pair_times(mask_list)
    )
}

#' Check that every mask_list element is a (subject, time) pair
#'
#' @description
#' A length-one element used to read `x[2]` as `NA` and fail much later with
#' an unrelated message, so the shape is checked up front. Note that
#' `list(A1 = 3)` is not such a pair: the subject belongs inside the element,
#' as `list(c("A1", 3))`.
#'
#' @param mask_list The list supplied by the user.
#'
#' @return `NULL`, invisibly. Called for the error it raises.
#'
#' @keywords internal
#' @noRd
mc_check_pair_shape <- function(mask_list) {
    len <- vapply(mask_list, length, integer(1))
    if (length(len) == 0 || !any(len != 2)) {
        return(invisible(NULL))
    }

    bad <- which(len != 2)
    stop(
        "Every element of mask_list must be a (subject, time) pair of ",
        "length 2. Element(s) ", mc_fmt_some(bad), " have length ",
        mc_fmt_some(unique(len[bad])), ". Write the pairs as ",
        "list(c(\"A1\", 3), c(\"A2\", 5)) -- note that ",
        "list(A1 = 3) is a named number, not a pair -- or supply a ",
        "data frame such as data.frame(rep = \"A1\", time = 3).",
        call. = FALSE
    )
}

#' Extract the time from each mask_list pair
#'
#' @param mask_list The list supplied by the user, already shape-checked.
#'
#' @return Numeric vector of times, one per element.
#'
#' @keywords internal
#' @noRd
mc_pair_times <- function(mask_list) {
    raw <- vapply(mask_list, function(x) as.character(x[2]), character(1))
    num <- as.numeric(raw[!is.na(raw) & mc_looks_numeric(raw)])

    if (length(num) != length(raw)) {
        bad <- which(is.na(raw) | !mc_looks_numeric(raw))
        stop(
            "The time in mask_list element(s) ", mc_fmt_some(bad),
            " is not numeric: ", mc_fmt_some(raw[bad]), ".",
            call. = FALSE
        )
    }
    num
}

#' Test whether strings would parse as numbers
#'
#' @description
#' Checked with a pattern rather than by converting, so that an invalid entry
#' reports itself instead of surfacing as a coercion warning.
#'
#' @param x Character vector to test.
#'
#' @return Logical vector, `TRUE` where `x` would parse as a number.
#'
#' @keywords internal
#' @noRd
mc_looks_numeric <- function(x) {
    pat <- "^[+-]?(([0-9]+[.]?[0-9]*)|([.][0-9]+))([eE][+-]?[0-9]+)?$"
    !is.na(x) & grepl(pat, trimws(x))
}

#' Turn an observation matrix into subject-timepoint pairs
#'
#' @param mask_matrix Matrix or data frame of 0/1 flags, subjects in rows and
#'   time points in columns. A `subject` column is used as row names when
#'   present.
#'
#' @return Data frame with `rep` and `time` columns, one row per zero flag.
#'
#' @keywords internal
#' @noRd
mc_mask_pairs_from_matrix <- function(mask_matrix) {
    mm <- as.data.frame(mask_matrix)

    if ("subject" %in% names(mm)) {
        subj <- mm$subject
        mm <- mm[, setdiff(names(mm), "subject"), drop = FALSE]
        rownames(mm) <- subj
    }

    mm_long <- tibble::rownames_to_column(mm, "rep")
    mm_long <- tidyr::pivot_longer(
        mm_long, -rep,
        names_to = "time", values_to = "obs_flag"
    )
    mm_long <- dplyr::mutate(mm_long, time = as.numeric(time))

    mm_missing <- dplyr::filter(mm_long, obs_flag == 0)
    dplyr::select(mm_missing, rep, time)
}

#' Validate a table and return its parsed sample columns
#'
#' @param dat Wide abundance table.
#' @param taxon_col Character name of the taxon column.
#' @param parse_fun Function mapping column names to subject and time.
#'
#' @return Data frame with `col`, `rep` and `time`, holding only the columns
#'   that parsed.
#'
#' @keywords internal
#' @noRd
mc_prepare_col_map <- function(dat, taxon_col, parse_fun) {
    if (!is.data.frame(dat)) {
        stop("dat must be a data.frame")
    }
    if (!(taxon_col %in% colnames(dat))) {
        stop("taxon_col not found in dat")
    }

    col_map <- parse_fun(setdiff(colnames(dat), taxon_col))
    if (!all(c("col", "rep", "time") %in% names(col_map))) {
        stop("parse_fun must return columns: col, rep, time")
    }

    col_map <- dplyr::filter(col_map, !is.na(rep), !is.na(time))
    if (nrow(col_map) == 0) {
        stop("No valid replicate/time columns detected.")
    }

    col_map
}
