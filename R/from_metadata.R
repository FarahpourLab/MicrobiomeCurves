# Accepting data in the shape a study actually produces.
#
# The package works internally on a wide table whose columns encode subject
# and time in their names. That is convenient for the fitting code and
# unreasonable to ask of anyone preparing data, who will have an abundance
# matrix (taxa in rows, samples in columns) and a separate metadata table
# saying which sample belongs to which subject at which time.
#
# This file translates between the two, reports what the design is missing,
# and keeps the mapping so results can be handed back under the sample names
# the user supplied. The encoding itself is shared with the
# SummarizedExperiment path, so both entry points agree.

#' Build a design from an abundance table and its metadata
#'
#' @description
#' Takes data in the form a study produces — an abundance table with taxa in
#' rows and samples in columns, plus a metadata table saying which sample
#' belongs to which subject at which time — and prepares it for imputation.
#'
#' The three metadata columns are named by argument, so they can be called
#' anything. `subject` refers to whatever the repeated measurements are taken
#' on: a mouse, a participant, a plot, a bioreactor.
#'
#' What the design is missing is reported as it is found, and returned in the
#' result. Two kinds of gap are distinguished: a sample listed in the
#' metadata whose abundance column holds no data, and a subject-timepoint
#' with no sample at all, which is discovered from the grid of subjects
#' against time points.
#'
#' @param abundance A data.frame or matrix with taxa in rows and samples in
#'   columns. Sample names are taken from the column names.
#' @param metadata A data.frame with one row per sample.
#' @param sample_col Character name of the metadata column holding sample
#'   identifiers. These must match the columns of `abundance`.
#' @param subject_col Character name of the metadata column identifying the
#'   subject each sample was taken from.
#' @param time_col Character name of the metadata column holding the time
#'   point. Values must be numeric or coercible to numeric.
#' @param taxon_col Character name of the column of `abundance` holding taxon
#'   identifiers, or `NULL` (the default) to take them from the row names.
#' @param verbose Logical. Whether to report the design as it is built.
#'
#' @return An object of class `tti_design`: a list with the converted
#'   `table`, the `map` from internal column names back to sample names, the
#'   `metadata` extended with any subject-timepoint that has no sample,
#'   `missing`, `observed`, `subjects`, `times` and `taxon_col`.
#'
#' @examples
#' counts <- data.frame(
#'     taxon = c("Bacteroides", "Prevotella"),
#'     S1 = c(1.2, 0.4), S2 = c(1.5, 0.6), S3 = c(1.1, 0.5)
#' )
#' meta <- data.frame(
#'     sample = c("S1", "S2", "S3"),
#'     mouse = c("M01", "M01", "M02"),
#'     day = c(0, 7, 0)
#' )
#'
#' design <- tti_from_metadata(
#'     counts, meta,
#'     sample_col = "sample", subject_col = "mouse", time_col = "day",
#'     taxon_col = "taxon"
#' )
#' design$missing
#'
#' @seealso [tti_run()], which accepts the same arguments and returns results
#'   under the original sample names.
#'
#' @export
tti_from_metadata <- function(abundance, metadata,
                              sample_col, subject_col, time_col,
                              taxon_col = NULL, verbose = TRUE) {
    say <- function(...) if (isTRUE(verbose)) message(...)

    ab <- tti_abundance_parts(abundance, taxon_col)
    md <- tti_metadata_parts(
        metadata, sample_col, subject_col, time_col, colnames(ab$mat)
    )

    design <- tti_build_design(ab, md, say)
    tti_report_design(design, say)
    design
}
