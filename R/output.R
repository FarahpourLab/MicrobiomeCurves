# Writing a run to disk.
#
# What is written is the record of the run: the design, the time points in
# order, every subject-timepoint that was missing and why, any warning
# raised, and the completed table itself. Someone opening the folder months
# later should be able to see what was imputed and on what basis, without
# rerunning anything.

#' Write a run's log and completed table to a directory
#'
#' @param run The `mc_run` object.
#' @param out_dir Directory to write into. Created if it does not exist.
#' @param warnings Character vector of warnings raised during the run.
#'
#' @return Character vector of the two paths written, invisibly.
#'
#' @keywords internal
#' @noRd
mc_write_output <- function(run, out_dir, warnings = character(0)) {
    if (!dir.exists(out_dir)) {
        ok <- dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        if (!ok) {
            stop(
                "Could not create out_dir '", out_dir, "'.",
                call. = FALSE
            )
        }
    }

    written <- mc_write_tables(run, out_dir)

    log <- file.path(out_dir, "imputation_log.txt")
    writeLines(mc_log_lines(run, warnings), log)

    invisible(c(written, log))
}

#' Write the completed table on every scale
#'
#' @description
#' The model works in centred log ratios, so that table is always written.
#' Relative abundance and counts are derived from it, because most
#' downstream analysis wants one of those.
#'
#' @param run The `mc_run` object.
#' @param out_dir Directory to write into.
#'
#' @return Character vector of the paths written.
#'
#' @keywords internal
#' @noRd
mc_write_tables <- function(run, out_dir) {
    info <- run$scale_info
    if (is.null(info)) info <- list(depths = NULL, pseudocount = NULL)

    ra <- mc_completed_ra(run$completed, mc_fill_scale(
        info$pseudocount, info$default_pseudocount, names(run$completed)
    ))
    known <- info$observed
    if (!is.null(known) && ncol(known) > 0) {
        ra <- mc_use_known(ra, sweep(known, 2, colSums(known), "/"))
    }
    depth <- info$default_depth
    if (is.null(known) || ncol(known) == 0) {
        # No library sizes were supplied. One depth is implied for the whole
        # table rather than one per sample, so that counts stay comparable
        # between samples.
        depth <- mc_implied_depth(as.matrix(ra[, -1, drop = FALSE]))
    }
    counts <- mc_completed_counts(ra, info$depths, depth)
    counts <- mc_use_known(counts, known)

    tables <- list(
        imputed_clr = run$completed,
        imputed_relative_abundance = ra,
        imputed_counts = counts
    )

    paths <- character(0)
    for (nm in names(tables)) {
        f <- file.path(out_dir, paste0(nm, ".tsv"))
        # Tab separated and unquoted. Taxonomy strings carry commas often
        # and tabs never, so the separator does the work a quote would.
        utils::write.table(
            tables[[nm]], f,
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
        )
        paths <- c(paths, f)
    }
    paths
}

#' Extend a per-sample value to every column of the completed table
#'
#' @description
#' A sample that was created has no value of its own, so it takes the
#' typical one. Without this it would fall through to a value derived from
#' its own numbers and land on a different scale from the study.
#'
#' @param values Named numeric vector over the caller's samples, or `NULL`.
#' @param default Value for a sample that has none.
#' @param cols Column names of the completed table.
#'
#' @return A named numeric vector over every sample column, or `NULL`.
#'
#' @keywords internal
#' @noRd
mc_fill_scale <- function(values, default, cols) {
    if (is.null(values)) {
        return(NULL)
    }
    samples <- setdiff(cols, "taxon")
    out <- stats::setNames(rep(default, length(samples)), samples)
    shared <- intersect(names(values), samples)
    out[shared] <- values[shared]
    out[is.na(out)] <- default
    out
}

#' Compose the lines of the run log
#'
#' @param run The `mc_run` object.
#' @param warnings Character vector of warnings raised during the run.
#'
#' @return Character vector of lines.
#'
#' @keywords internal
#' @noRd
mc_log_lines <- function(run, warnings) {
    d <- run$design

    c(
        "MicrobiomeCurves run log",
        paste0("package version: ", mc_version_string()),
        "",
        "DESIGN",
        paste0("  taxa        : ", nrow(run$completed)),
        paste0("  samples     : ", nrow(d$map)),
        paste0("  subjects    : ", length(d$subjects)),
        paste0("  time points : ", length(d$times)),
        paste0("  time order  : ", mc_axis_text(d$axis)),
        paste0("  subjects    : ", mc_fmt_some(d$subjects, n = 20)),
        "",
        mc_log_missing(run, d),
        "",
        mc_log_result(run),
        "",
        mc_log_scales(run),
        "",
        mc_log_warnings(warnings)
    )
}

#' The missing-sample section of the log
#'
#' @param run The `mc_run` object.
#' @param d Its design.
#'
#' @return Character vector of lines.
#'
#' @keywords internal
#' @noRd
mc_log_missing <- function(run, d) {
    n_cells <- length(d$subjects) * length(d$times)
    miss <- d$missing

    head <- paste0(
        "MISSING (", nrow(miss), " of ", n_cells, " subject-timepoints)"
    )
    if (nrow(miss) == 0) {
        return(c(head, "  none: every subject has a sample at every time."))
    }

    body <- vapply(
        seq_len(nrow(miss)),
        function(i) {
            paste0(
                "  ", miss$subject[i], "  at time ", miss$time_label[i],
                "  [", miss$reason[i], "]",
                if (!is.na(miss$sample[i])) {
                    paste0("  sample ", miss$sample[i])
                } else {
                    "  no sample was collected"
                }
            )
        },
        character(1)
    )
    c(head, body)
}

#' The outcome section of the log
#'
#' @param run The `mc_run` object.
#'
#' @return Character vector of lines.
#'
#' @keywords internal
#' @noRd
mc_log_result <- function(run) {
    n_cells <- nrow(run$imputed)
    created <- sum(run$metadata$imputed)

    c(
        "RESULT",
        paste0("  cells imputed  : ", n_cells - run$n_failed, " of ", n_cells),
        paste0("  columns added  : ", created),
        paste0(
            "  added columns  : ",
            if (created > 0) {
                mc_fmt_some(run$metadata$sample[run$metadata$imputed], n = 20)
            } else {
                "none"
            }
        )
    )
}

#' The output-scale section of the log
#'
#' @param run The `mc_run` object.
#'
#' @return Character vector of lines.
#'
#' @keywords internal
#' @noRd
mc_log_scales <- function(run) {
    raw <- !is.null(run$scale_info$observed)

    c(
        "OUTPUT SCALES",
        "  imputed_clr.tsv                : what the model produced",
        "  imputed_relative_abundance.tsv : proportions, summing to one",
        "  imputed_counts.tsv             : proportions times a library size",
        "",
        if (raw) {
            c(
                "  You supplied counts, so a sample that was collected is",
                "  written back exactly as you gave it. A sample that was",
                "  created carries the median library size of the study, and",
                "  its zeros are the taxa falling below the value the CLR",
                "  step used in place of a zero."
            )
        } else {
            c(
                "  You supplied CLR values, so library sizes and structural",
                "  zeros were not available. Proportions come from inverting",
                "  the CLR and are strictly positive: a zero in your original",
                "  table appears here as a small positive number. Counts use",
                "  a library size implied by the rarest taxon in each sample,",
                "  so they are on a plausible scale, not the original one."
            )
        }
    )
}

#' The warnings section of the log
#'
#' @param warnings Character vector of warnings raised during the run.
#'
#' @return Character vector of lines.
#'
#' @keywords internal
#' @noRd
mc_log_warnings <- function(warnings) {
    if (length(warnings) == 0) {
        return(c("WARNINGS", "  none"))
    }
    c(
        paste0("WARNINGS (", length(warnings), ")"),
        paste0("  - ", gsub("[[:space:]]+", " ", warnings))
    )
}

#' The package version, for the log header
#'
#' @return A single string.
#'
#' @keywords internal
#' @noRd
mc_version_string <- function() {
    v <- tryCatch(
        as.character(utils::packageVersion("MicrobiomeCurves")),
        error = function(e) "unknown"
    )
    v
}
