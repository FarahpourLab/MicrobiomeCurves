# Running one study as several independent studies.
#
# A trial with treatment arms, or a cohort split by disease status, should
# not borrow a curve shape from one group to fill a gap in another. When the
# caller names a group column, the samples are split on it and each
# group is imputed on its own: its own outlier screening, its own
# clustering, its own FPCA. Nothing crosses between groups.
#
# Every group writes the full set of files it would have written as a run of
# its own, into a subdirectory named after it, and the completed tables are
# also merged back into the caller's original layout at the top level.

#' Group sample names by the group their subject belongs to
#'
#' @description
#' Splits the samples on a group column, after checking that the split
#' is well defined.
#'
#' @details
#' A sample that is missing entirely has no row in the metadata, so there is
#' nothing to read a group from. Its group can only be inherited
#' from its subject, which means a subject has to sit in exactly one
#' group for the whole study. A subject whose group changes between
#' time points would be cut in half by the split, leaving too few
#' observations to fit a curve to either piece, so that is an error rather
#' than a warning.
#'
#' @param metadata The metadata data frame.
#' @param sample_col,subject_col,group_col Column names within it.
#'
#' @return Named list of character vectors of sample names, one entry per
#'   group, in order of first appearance.
#'
#' @keywords internal
#' @noRd
mc_group_samples <- function(metadata, sample_col, subject_col,
                                 group_col) {
    if (!is.character(group_col) || length(group_col) != 1 ||
        is.na(group_col)) {
        stop("'group_col' must be a single column name.", call. = FALSE)
    }
    if (!group_col %in% names(metadata)) {
        stop(
            "Column '", group_col, "' is not in the metadata. ",
            "Available: ", paste(names(metadata), collapse = ", "), ".",
            call. = FALSE
        )
    }

    cond <- as.character(metadata[[group_col]])
    subj <- as.character(metadata[[subject_col]])

    if (anyNA(cond)) {
        bad <- unique(subj[is.na(cond)])
        stop(
            "Column '", group_col, "' is missing for ", sum(is.na(cond)),
            " sample(s), in subject(s): ", mc_short_list(bad), ". ",
            "Every sample needs a group before the run can be split.",
            call. = FALSE
        )
    }

    # One subject, one group. See the note above on why this cannot be
    # relaxed into a warning.
    per_subject <- tapply(cond, subj, function(x) length(unique(x)))
    mixed <- names(per_subject)[per_subject > 1]
    if (length(mixed) > 0) {
        stop(
            "Subject(s) ", mc_short_list(mixed), " appear under more than ",
            "one value of '", group_col, "'. A group has to hold ",
            "for a subject's whole time course, because a sample that is ",
            "missing altogether has no metadata row to read one from.",
            call. = FALSE
        )
    }

    groups <- split(as.character(metadata[[sample_col]]), factor(cond, unique(cond)))
    lapply(groups, as.character)
}

#' Format a character vector for an error message
#'
#' @param x Character vector.
#' @param n Most to show before abbreviating.
#'
#' @return A single string.
#'
#' @keywords internal
#' @noRd
mc_short_list <- function(x, n = 5) {
    x <- unique(x)
    if (length(x) <= n) {
        return(paste(x, collapse = ", "))
    }
    paste0(paste(utils::head(x, n), collapse = ", "),
           " and ", length(x) - n, " more")
}

#' Take the columns of one group out of an abundance table
#'
#' @param dat The abundance table, a matrix or a data frame.
#' @param keep Character vector of sample names to retain.
#'
#' @return `dat`, holding only those samples, in their original order.
#'
#' @keywords internal
#' @noRd
mc_subset_samples <- function(dat, keep) {
    cols <- if (is.matrix(dat)) colnames(dat) else names(dat)
    sel <- cols[cols %in% keep]
    if (is.matrix(dat)) dat[, sel, drop = FALSE] else dat[, sel, drop = FALSE]
}

#' Impute each group separately and bring the results back together
#'
#' @description
#' The body of [mc_run()] when `group_col` is given.
#'
#' @details
#' Each group is a complete run in its own right, written to its own
#' subdirectory. A group that cannot be fitted, or has nothing missing,
#' does not stop the others: it is reported and its samples are returned as
#' they came in.
#'
#' Scale conversion is the one thing left global. Pseudocount and library
#' depth describe the table that was handed in, not the model, and deriving
#' them per group would put the count tables of two groups on
#' different scales and make them incomparable.
#'
#' @param dat,metadata,sample_col,subject_col,time_col As in [mc_run()].
#' @param abundance_type,pseudocount,out_dir As in
#'   [mc_run()].
#' @param C,cluster_method,use_outliers,seed,min_observed,verbose As in
#'   [mc_run()].
#' @param group_col Name of the group column in `metadata`.
#'
#' @return An object of class `mc_run`, with `groups` holding the
#'   individual runs.
#'
#' @keywords internal
#' @noRd
mc_run_groups <- function(dat, metadata, sample_col, subject_col,
                               time_col, abundance_type, pseudocount,
                               out_dir,
                               C, cluster_method, use_outliers,
                               seed, min_observed, verbose, group_col) {
    groups <- mc_group_samples(
        metadata, sample_col, subject_col, group_col
    )
    say <- function(...) if (isTRUE(verbose)) message(...)
    say(
        "Splitting on '", group_col, "': ", length(groups),
        " group(s) run separately (", paste(names(groups), collapse = ", "),
        ")."
    )

    runs <- list()
    notes <- character(0)
    for (nm in names(groups)) {
        say("--- group '", nm, "' ---")
        sub_dir <- if (is.null(out_dir)) NULL else {
            file.path(out_dir, mc_safe_name(nm))
        }
        res <- mc_try_group(
            dat, metadata, groups[[nm]], nm, sample_col, subject_col,
            time_col, abundance_type, pseudocount, sub_dir,
            C, cluster_method, use_outliers, seed,
            min_observed, verbose
        )
        if (is.null(res$run)) {
            notes <- c(notes, res$note)
            say("    skipped: ", res$note)
        } else {
            runs[[nm]] <- res$run
        }
    }

    if (length(runs) == 0) {
        stop(
            "No group could be imputed. ",
            paste(notes, collapse = " "),
            call. = FALSE
        )
    }

    merged <- mc_merge_groups(runs, dat, group_col)
    merged <- mc_attach_scale_info(merged, dat, abundance_type)
    merged$groups <- runs
    merged$group_col <- group_col

    warned <- c(notes, unlist(lapply(runs, function(r) r$warnings)))
    merged$warnings <- unique(as.character(warned))

    if (!is.null(out_dir)) {
        # The merged tables only. Each group has already written its own
        # log, which cannot be merged because each describes a different
        # set of subjects.
        written <- mc_write_output(merged, out_dir, merged$warnings)
        merged$files <- c(
            written, unlist(lapply(runs, function(r) r$files), use.names = FALSE)
        )
        if (isTRUE(verbose)) {
            message("Wrote merged tables to ", out_dir, ".")
        }
    }

    merged
}

#' Run one group, turning a failure into a note
#'
#' @param dat,metadata The full inputs.
#' @param keep Character vector of this group's sample names.
#' @param nm The group's label.
#' @param sample_col,subject_col,time_col,abundance_type,pseudocount As in
#'   [mc_run()].
#' @param sub_dir Directory for this group, or `NULL`.
#' @param C,cluster_method,use_outliers,seed As in
#'   [mc_run()].
#' @param min_observed,verbose As in [mc_run()].
#'
#' @return List with `run` (or `NULL`) and `note`.
#'
#' @keywords internal
#' @noRd
mc_try_group <- function(dat, metadata, keep, nm, sample_col,
                              subject_col, time_col, abundance_type,
                              pseudocount, sub_dir,
                              C, cluster_method, use_outliers, seed,
                              min_observed, verbose) {
    md <- metadata[as.character(metadata[[sample_col]]) %in% keep, ,
                   drop = FALSE]
    ab <- mc_subset_samples(dat, keep)

    run <- tryCatch(
        mc_run_from_metadata(
            ab, md, sample_col, subject_col, time_col, abundance_type,
            pseudocount, sub_dir, C,
            cluster_method, use_outliers, seed, min_observed, verbose
        ),
        error = function(e) {
            structure(conditionMessage(e), class = "mc_group_error")
        }
    )

    if (inherits(run, "mc_group_error")) {
        return(list(
            run = NULL,
            note = paste0(
                "Group '", nm, "' was not imputed: ", as.character(run)
            )
        ))
    }
    list(run = run, note = NULL)
}

#' Put the groups' completed tables back into one table
#'
#' @description
#' Each group filled its own samples. This writes all of them back into
#' the layout the caller handed in, so the merged table has the original
#' columns in the original order, plus a column for every sample that was
#' created because it had been missing entirely.
#'
#' @param runs Named list of finished runs, one per group.
#' @param dat The abundance table as supplied.
#' @param group_col Name of the group column, for labelling.
#'
#' @return An `mc_run` holding the merged tables.
#'
#' @keywords internal
#' @noRd
mc_merge_groups <- function(runs, dat, group_col) {
    taxa <- mc_merge_taxa(runs)

    # Start from the table as supplied. A group that could not be
    # imputed still has its samples here, carried through untouched, rather
    # than disappearing from the merged output.
    completed <- mc_original_table(dat, taxa)

    for (nm in names(runs)) {
        cmp <- runs[[nm]]$completed
        idx <- match(taxa, cmp[[1]])
        for (col in setdiff(names(cmp), names(cmp)[1])) {
            completed[[col]] <- cmp[[col]][idx]
        }
    }

    # Original column order first, then whatever each group created.
    original <- if (is.matrix(dat)) colnames(dat) else names(dat)
    original <- intersect(original, names(completed))
    rest <- setdiff(names(completed), c("taxon", original))
    completed <- completed[, c("taxon", original, rest), drop = FALSE]

    imputed <- mc_merge_imputed(runs, group_col)
    design <- mc_merge_design(runs)

    structure(
        list(
            completed = completed,
            imputed = imputed,
            taxon_col = "taxon",
            n_failed = sum(vapply(
                runs, function(r) as.integer(r$n_failed %||% 0L), integer(1)
            )),
            design = design,
            metadata = do.call(rbind, lapply(runs, function(r) r$metadata)),
            missing = design$missing,
            observed = NULL
        ),
        class = "mc_run"
    )
}

#' Combine the groups' designs into one, for the merged log
#'
#' @description
#' Subjects are disjoint across groups, by the one-subject-one-group
#' rule, so concatenating them is unambiguous. Time points are shared, so
#' they are taken in their first group's order and extended only if a
#' later one adds to them.
#'
#' @param runs Named list of finished runs.
#'
#' @return A list shaped like an `mc_design`, holding only what the log
#'   reads.
#'
#' @keywords internal
#' @noRd
mc_merge_design <- function(runs) {
    first <- runs[[1]]$design

    times <- first$times
    for (r in runs) times <- union(times, r$design$times)
    if (is.numeric(first$times)) times <- sort(as.numeric(times))

    list(
        map = do.call(rbind, lapply(runs, function(r) r$design$map)),
        subjects = unlist(lapply(runs, function(r) r$design$subjects),
                          use.names = FALSE),
        times = times,
        axis = first$axis,
        taxon_col = first$taxon_col,
        missing = do.call(rbind, lapply(runs, function(r) r$design$missing))
    )
}

#' The table as supplied, as a taxon column plus its sample columns
#'
#' @param dat The abundance table, a matrix or a data frame.
#' @param taxa Character vector of taxon names to carry, in output order.
#'
#' @return A data frame with a `taxon` column and one column per sample.
#'
#' @keywords internal
#' @noRd
mc_original_table <- function(dat, taxa) {
    out <- data.frame(taxon = taxa, stringsAsFactors = FALSE)

    if (is.matrix(dat)) {
        idx <- match(taxa, rownames(dat))
        for (col in colnames(dat)) out[[col]] <- dat[idx, col]
        return(out)
    }

    idx <- match(taxa, as.character(dat[[1]]))
    for (col in names(dat)[-1]) out[[col]] <- dat[[col]][idx]
    out
}

#' The taxa of every group, in the order the first one had them
#'
#' @param runs Named list of finished runs.
#'
#' @return Character vector of taxon names.
#'
#' @keywords internal
#' @noRd
mc_merge_taxa <- function(runs) {
    out <- NULL
    for (r in runs) {
        this <- as.character(r$completed[[1]])
        out <- if (is.null(out)) this else union(out, this)
    }
    out
}

#' Stack the groups' long tables, labelled by group
#'
#' @param runs Named list of finished runs.
#' @param group_col Name to give the label column.
#'
#' @return One data frame, or `NULL` when there is nothing to stack.
#'
#' @keywords internal
#' @noRd
mc_merge_imputed <- function(runs, group_col) {
    parts <- list()
    for (nm in names(runs)) {
        p <- runs[[nm]]$imputed
        if (is.null(p) || nrow(p) == 0) next
        p <- as.data.frame(p, stringsAsFactors = FALSE)
        p[[group_col]] <- nm
        parts[[nm]] <- p
    }
    if (length(parts) == 0) {
        return(NULL)
    }

    keep <- Reduce(intersect, lapply(parts, names))
    do.call(rbind, lapply(parts, function(p) p[, keep, drop = FALSE]))
}

`%||%` <- function(a, b) if (is.null(a)) b else a
