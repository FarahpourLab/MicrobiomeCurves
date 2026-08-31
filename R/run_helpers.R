#' Check the table layout and return its usable sample columns
#'
#' @param dat Wide abundance table.
#' @param taxon_col Character name of the taxon column.
#' @param parse_fun Function mapping column names to subject and time.
#'
#' @return List with `col_map`, the parsed sample columns, and `unparsed`,
#'   the column names that did not match.
#'
#' @keywords internal
#' @noRd
mc_column_map <- function(dat, taxon_col, parse_fun) {
    if (!is.data.frame(dat)) {
        stop("dat must be a data.frame")
    }
    if (!(taxon_col %in% colnames(dat))) {
        stop("taxon_col '", taxon_col, "' not found in dat")
    }

    data_cols <- setdiff(colnames(dat), taxon_col)
    col_map <- parse_fun(data_cols)

    if (!all(c("col", "rep", "time") %in% names(col_map))) {
        stop("parse_fun must return columns: col, rep, time")
    }

    unparsed <- col_map$col[is.na(col_map$rep) | is.na(col_map$time)]
    keep_cols <- !is.na(col_map$rep) & !is.na(col_map$time)
    col_map <- col_map[keep_cols, , drop = FALSE]

    if (nrow(col_map) == 0) {
        stop(
            "No '<subject>.<time>' columns detected. Expected names such as ",
            "'B005.0'. Supply parse_fun if your table uses another convention."
        )
    }

    # A column that is entirely NA carries no data, so its storage type does
    # not matter: read.csv() returns such a column as logical, and that is a
    # normal way for a missing sample to arrive.
    usable <- vapply(
        dat[col_map$col],
        function(x) is.numeric(x) || all(is.na(x)),
        logical(1)
    )
    non_numeric <- col_map$col[!usable]

    if (length(non_numeric) > 0) {
        stop(
            "These sample columns are not numeric: ",
            paste(utils::head(non_numeric, 5), collapse = ", "),
            if (length(non_numeric) > 5) ", ..." else "",
            "\nAbundances must be numeric",
            " (the method expects CLR-transformed values)."
        )
    }

    list(col_map = col_map, unparsed = unparsed)
}

#' Find the samples that are missing from a table
#'
#' @param dat Wide abundance table.
#' @param col_map Parsed sample columns.
#' @param reps Character vector of subjects.
#' @param times Numeric vector of the intended time grid.
#' @param make_col Function of `(rep, time)` returning a column name.
#'
#' @return List with `missing`, the missing subject-timepoints and why, and
#'   `partial_na`, the columns that hold some but not only `NA`.
#'
#' @keywords internal
#' @noRd
mc_find_missing <- function(dat, col_map, reps, times, make_col) {
    n_na <- vapply(dat[col_map$col], function(x) sum(is.na(x)), numeric(1))
    n_taxa <- nrow(dat)
    all_na <- col_map$col[n_na == n_taxa]

    partial_na <- mc_partial_na_table(dat, col_map)

    missing_all_na <- as.data.frame(
        col_map[col_map$col %in% all_na, c("rep", "time"), drop = FALSE],
        stringsAsFactors = FALSE
    )
    if (nrow(missing_all_na) > 0) missing_all_na$reason <- "all_na"

    grid <- expand.grid(
        rep = reps, time = times,
        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
    present_key <- paste(col_map$rep, col_map$time, sep = "\r")
    grid_key <- paste(grid$rep, grid$time, sep = "\r")

    missing_absent <- grid[!(grid_key %in% present_key), , drop = FALSE]
    if (nrow(missing_absent) > 0) missing_absent$reason <- "absent_column"

    missing <- rbind(missing_absent, missing_all_na)

    if (nrow(missing) > 0) {
        # radix ordering uses C collation, so row order is identical on every
        # platform and locale. Default collation is locale-dependent.
        ord <- order(missing$rep, missing$time, method = "radix")
        missing <- missing[ord, , drop = FALSE]
        rownames(missing) <- NULL
        missing$col <- make_col(missing$rep, missing$time)
    } else {
        missing <- data.frame(
            rep = character(0), time = numeric(0),
            reason = character(0), col = character(0),
            stringsAsFactors = FALSE
        )
    }

    list(missing = missing, partial_na = partial_na)
}

#' Tabulate the columns that are only partly measured
#'
#' @description
#' A taxon that is `NA` everywhere accounts for its own missing values, so its
#' row is excluded before deciding whether a column is only partly measured.
#' Otherwise a single empty taxon would flag every sample in the table.
#'
#' @param dat The user's data frame.
#' @param col_map Data frame of parsed columns, with a `col` column.
#'
#' @return Data frame with `col`, `n_na` and `n_taxa`, with no rows when every
#'   column is either fully observed or entirely `NA`.
#'
#' @keywords internal
#' @noRd
mc_partial_na_table <- function(dat, col_map) {
    keep <- mc_informative_taxa(dat, col_map)
    n_keep <- sum(keep)

    n_na <- if (n_keep == 0) {
        stats::setNames(rep(0, length(col_map$col)), col_map$col)
    } else {
        vapply(
            dat[keep, col_map$col, drop = FALSE],
            function(x) sum(is.na(x)), numeric(1)
        )
    }

    partial <- col_map$col[n_na > 0 & n_na < n_keep]
    data.frame(
        col = partial,
        n_na = as.numeric(n_na[partial]),
        n_taxa = rep(n_keep, length(partial)),
        stringsAsFactors = FALSE
    )
}

#' Count how many time points each subject actually has
#'
#' @param missing Data frame of absent (subject, time) pairs.
#' @param reps Character vector of subject identifiers.
#' @param times Numeric vector of the time points in the design.
#'
#' @return Data frame with `rep` and `n_observed`.
#'
#' @keywords internal
#' @noRd
mc_observed_counts <- function(missing, reps, times) {
    missing_key <- paste(missing$rep, missing$time, sep = "\r")
    observed_n <- vapply(
        reps,
        function(r) sum(!(paste(r, times, sep = "\r") %in% missing_key)),
        numeric(1)
    )

    data.frame(
        rep = reps,
        n_observed = as.numeric(observed_n),
        stringsAsFactors = FALSE
    )
}

#' Assemble the object returned by mc_run()
#'
#' @param res List from [mc_run_pipeline()].
#' @param info List from [mc_survey_input()].
#' @param taxon_col Character name of the taxon column.
#'
#' @return An object of class `mc_run`.
#'
#' @keywords internal
#' @noRd
mc_run_result <- function(res, info, taxon_col) {
    structure(
        list(
            completed = res$completed,
            imputed = res$pred,
            missing = info$missing,
            observed = info$observed,
            partial_na = info$partial_na,
            n_failed = res$n_failed,
            taxon_col = taxon_col,
            fit = res$fit
        ),
        class = "mc_run"
    )
}

#' Warn about anything in the table that needs the user's attention
#'
#' @param info List returned by [mc_detect_missing()].
#' @param min_observed Integer. Observation count below which a subject is
#'   flagged.
#'
#' @return `NULL`, invisibly. Called for its warnings.
#'
#' @keywords internal
#' @noRd
mc_warn_about_input <- function(info, min_observed,
                                subject_label = identity) {
    # The column warnings are raised by mc_detect_missing(), which this path
    # has already been through, so only the thin-subject check runs here.
    thin <- info$observed$rep[info$observed$n_observed < min_observed]
    if (length(thin) > 0) {
        # The fit works in internal subject codes. Name the subjects the way
        # the caller wrote them, or the warning is unactionable.
        named <- subject_label(thin)
        warning(
            length(named), " subject(s) have fewer than ", min_observed,
            " observed time points: ",
            paste(utils::head(named, 5), collapse = ", "),
            if (length(named) > 5) ", ..." else "",
            ". Their imputations fall back towards the population mean ",
            "curve and carry little subject-specific information.",
            call. = FALSE
        )
    }

    invisible(NULL)
}

#' Detect what is missing and warn about anything unusual
#'
#' @param dat Wide abundance table.
#' @param taxon_col Character name of the taxon column.
#' @param times Optional numeric time grid.
#' @param parse_fun,make_col Column name parser and builder.
#' @param min_observed Integer threshold for flagging thin subjects.
#' @param say Function used to report progress.
#'
#' @return The list returned by [mc_detect_missing()].
#'
#' @keywords internal
#' @noRd
mc_survey_input <- function(dat, taxon_col, times, parse_fun, make_col,
                             min_observed, say, subject_label = identity) {
    info <- mc_detect_missing(
        dat = dat,
        taxon_col = taxon_col,
        times = times,
        parse_fun = parse_fun,
        make_col = make_col
    )

    if (nrow(info$missing) == 0) {
        stop(
            "No missing samples found: every subject-timepoint is ",
            "present and has at least one non-NA value. There is ",
            "nothing to impute."
        )
    }

    say(
        "Detected ", nrow(info$missing), " missing sample(s) across ",
        length(info$reps), " subject(s) and ",
        length(info$times), " time point(s)."
    )

    mc_warn_about_input(info, min_observed, subject_label)
    info
}

#' Give every missing cell a column to be written into
#'
#' @description
#' Creates an all-NA column for each subject-timepoint that has none, and
#' makes every sample column numeric.
#'
#' @param dat Wide abundance table.
#' @param info List returned by [mc_detect_missing()].
#' @param say Function used to report progress.
#'
#' @return The table, with columns added where needed.
#'
#' @keywords internal
#' @noRd
mc_add_missing_columns <- function(dat, info, say) {
    new_cols <- info$missing$col[info$missing$reason == "absent_column"]
    new_cols <- setdiff(unique(new_cols), colnames(dat))

    if (length(new_cols) > 0) {
        say(
            "Creating ", length(new_cols),
            " empty column(s) for absent samples."
        )
        for (cl in new_cols) {
            dat[[cl]] <- NA_real_
        }
    }

    # read.csv() gives an all-NA column logical type. Make every sample column
    # numeric so the table has one consistent type before reshaping.
    for (cl in info$col_map$col) {
        if (!is.numeric(dat[[cl]])) {
            dat[[cl]] <- as.numeric(dat[[cl]])
        }
    }

    dat
}

#' Write predictions back into the wide table
#'
#' @param dat Wide abundance table containing a column for every missing cell.
#' @param pred Long prediction table from [mc_impute()].
#'
#' @return The table with imputed values filled in.
#'
#' @keywords internal
#' @noRd
mc_fill_completed <- function(dat, pred) {
    for (cl in unique(pred$col)) {
        idx <- which(pred$col == cl)
        dat[[cl]][pred$taxon_idx[idx]] <- pred$imputed_value[idx]
    }
    dat
}

#' Report cells that could not be imputed
#'
#' @param pred Long prediction table from [mc_impute()].
#' @param say Function used to report progress.
#'
#' @return Integer count of cells left as `NA`.
#'
#' @keywords internal
#' @noRd
mc_report_failures <- function(pred, say) {
    # The warning itself is raised by mc_fit(), which this path has already
    # been through, so only the progress line is emitted here. Warning again
    # would report the same cells twice.
    n_failed <- sum(is.na(pred$imputed_value))

    say(
        "Done: ", nrow(pred) - n_failed, " of ", nrow(pred),
        " cell(s) imputed."
    )

    n_failed
}

#' Prepare, fit and write back the imputed values
#'
#' @param dat_full Table with a column for every missing cell.
#' @param taxon_col Character name of the taxon column.
#' @param parse_fun Column name parser.
#' @param info List returned by [mc_detect_missing()].
#' @param K Integer or `NULL`. Number of clusters.
#' @param cluster_method Character clustering method.
#' @param use_outliers Logical. Whether to screen outlying curves.
#' @param seed Integer random seed.
#' @param say Function used to report progress.
#'
#' @return List with the `fit`, the long predictions `pred`, the `completed`
#'   table and the `n_failed` count.
#'
#' @keywords internal
#' @noRd
mc_run_pipeline <- function(dat_full, taxon_col, parse_fun, info,
                             K, cluster_method, use_outliers, seed,
                             say) {
    prep <- mc_prepare(
        dat = dat_full,
        taxon_col = taxon_col,
        mask_list = info$missing[, c("rep", "time")],
        parse_fun = parse_fun
    )

    say("Fitting FPCA model over ", nrow(dat_full), " taxa ...")

    fit <- mc_fit(
        prep = prep,
        K = K,
        cluster_method = cluster_method,
        use_outliers = use_outliers,
        seed = seed
    )

    pred <- mc_impute(fit)

    list(
        fit = fit,
        pred = pred,
        completed = mc_fill_completed(dat_full, pred),
        n_failed = mc_report_failures(pred, say)
    )
}
