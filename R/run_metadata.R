# Running from an abundance table plus metadata.
#
# The caller never sees the internal column encoding. Their sample names go
# in, the design is reported, the fit runs on the encoded table, and the
# completed matrix comes back with their sample names restored. Columns
# created for a subject-timepoint that had no sample cannot carry a name of
# theirs, so one is built from the subject and the time and recorded in the
# returned metadata with imputed = TRUE.

#' Impute from an abundance table and its metadata
#'
#' @param dat The abundance table.
#' @param metadata The metadata data.frame.
#' @param sample_col,subject_col,time_col Column names within `metadata`.
#' @param taxon_col Character name of the taxon column of `dat`, or `NULL`
#'   to use row names.
#' @param K,cluster_method,use_outliers,seed Passed to the fit.
#' @param min_observed Integer. Subjects with fewer observed time points are
#'   reported.
#' @param verbose Logical. Whether to report progress.
#'
#' @return An object of class `tti_run`, with `completed` carrying the
#'   caller's sample names and an added `design` element.
#'
#' @keywords internal
#' @noRd
tti_run_from_metadata <- function(dat, metadata, sample_col, subject_col,
                                  time_col, taxon_col, K, cluster_method,
                                  use_outliers, seed, min_observed, verbose) {
    tti_check_meta_args(sample_col, subject_col, time_col)

    design <- tti_from_metadata(
        abundance = dat, metadata = metadata,
        sample_col = sample_col, subject_col = subject_col,
        time_col = time_col, taxon_col = taxon_col, verbose = verbose
    )

    run <- tti_run(
        design$table,
        taxon_col = design$taxon_col,
        K = K, cluster_method = cluster_method,
        use_outliers = use_outliers, seed = seed,
        min_observed = min_observed, verbose = verbose
    )

    tti_restore_names(run, design, taxon_col)
}

#' Require all three metadata column arguments together
#'
#' @param sample_col,subject_col,time_col The arguments as supplied.
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
#' @noRd
tti_check_meta_args <- function(sample_col, subject_col, time_col) {
    given <- c(
        sample_col = !is.null(sample_col),
        subject_col = !is.null(subject_col),
        time_col = !is.null(time_col)
    )
    if (all(given)) {
        return(invisible(NULL))
    }
    stop(
        "When metadata is supplied, sample_col, subject_col and time_col ",
        "must all be named. Missing: ",
        tti_fmt_some(names(given)[!given]), ".",
        call. = FALSE
    )
}

#' Put the caller's sample names back on a completed table
#'
#' @param run The `tti_run` object produced on the encoded table.
#' @param design The `tti_design` used to encode it.
#' @param taxon_col The caller's taxon column name, or `NULL`.
#'
#' @return The `tti_run` object, renamed, with `design` added.
#'
#' @keywords internal
#' @noRd
tti_restore_names <- function(run, design, taxon_col) {
    completed <- run$completed
    cols <- setdiff(names(completed), design$taxon_col)

    known <- match(cols, design$map$column)
    labels <- ifelse(
        is.na(known),
        tti_name_added(cols, design),
        design$map$sample[known]
    )

    id <- if (is.null(taxon_col)) "taxon" else taxon_col
    names(completed) <- c(id, labels)

    run$completed <- tti_order_samples(completed, cols, design, id)
    run$design <- design
    run$metadata <- tti_extend_metadata(design, cols[is.na(known)], labels)
    run$missing <- design$missing
    run$observed <- design$observed
    run
}

#' Put the completed table into subject and time order
#'
#' @description
#' A column created for an absent sample would otherwise land at the far
#' right, away from the rest of its subject. Ordering by subject and then
#' time puts every sample where a reader expects to find it.
#'
#' @param completed The completed table, already renamed.
#' @param cols Character vector of internal column names, in table order.
#' @param design The `tti_design` the run was built from.
#' @param id Name of the taxon identifier column.
#'
#' @return `completed`, with its sample columns reordered.
#'
#' @keywords internal
#' @noRd
tti_order_samples <- function(completed, cols, design, id) {
    parsed <- tti_parse_cols(cols)
    subject <- design$subjects[as.integer(sub("^s", "", parsed$rep))]
    time <- design$times[parsed$time + 1L]

    ord <- order(
        match(subject, design$subjects), time,
        method = "radix"
    )
    completed[, c(id, setdiff(names(completed), id)[ord]), drop = FALSE]
}
