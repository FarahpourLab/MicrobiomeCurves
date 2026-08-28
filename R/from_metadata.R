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
#' @param abundance A data.frame or matrix with taxa in the **row names** and
#'   one column per sample. Column names are the sample identifiers.
#' @param metadata A data.frame with one row per sample.
#' @param sample_col Character name of the metadata column holding sample
#'   identifiers. These must match the columns of `abundance`.
#' @param subject_col Character name of the metadata column identifying the
#'   subject each sample was taken from.
#' @param time_col Character name of the metadata column holding the time
#'   point. Either numbers or labels such as `"baseline"` and `"week1"`;
#'   labels take their order from the factor's levels if the column is a
#'   factor and from the row order otherwise, which is reported so it can be
#'   checked. **The model uses the order of the time points, not their
#'   values**, so uneven intervals are not represented in the curve; the
#'   values are kept for reporting and for the axis of the uncertainty
#'   pages.
#' @param abundance_type Either `"clr"`, meaning the values are already
#'   centred log-ratios, or `"raw"`, meaning counts or relative abundances to
#'   be CLR-transformed here.
#' @param pseudocount Used only when `abundance_type = "raw"`. Either
#'   `"auto"`, replacing zeros per sample with a fraction of that sample's
#'   smallest non-zero value, or a single positive number added to every
#'   entry.
#' @param verbose Logical. Whether to report the design as it is built.
#'
#' @return An object of class `mc_design`: a list with the converted
#'   `table`, the `map` from internal column names back to sample names, the
#'   `metadata`, `missing`, `observed`, `subjects`, `times`, the time `axis`
#'   and `taxon_col`.
#'
#' @examples
#' counts <- matrix(
#'     c(12, 4, 9, 15, 6, 11, 10, 5, 8),
#'     nrow = 3,
#'     dimnames = list(
#'         c("Bacteroides", "Faecalibacterium", "Bifidobacterium"),
#'         c("S001", "S002", "S003")
#'     )
#' )
#' meta <- data.frame(
#'     SampleID = c("S001", "S002", "S003"),
#'     SubjectID = c("SUB01", "SUB01", "SUB02"),
#'     Day = c(0, 7, 0)
#' )
#'
#' design <- mc_from_metadata(
#'     counts, meta,
#'     sample_col = "SampleID", subject_col = "SubjectID",
#'     time_col = "Day",
#'     abundance_type = "raw"
#' )
#' design$missing
#'
#' @seealso [mc_run()], which accepts the same arguments and returns results
#'   under the original sample names.
#'
#' @export
mc_from_metadata <- function(abundance, metadata,
                              sample_col, subject_col, time_col,
                              abundance_type = c("clr", "raw"),
                              pseudocount = "auto",
                              verbose = TRUE) {
    abundance_type <- match.arg(abundance_type)
    mc_check_pseudocount(pseudocount)
    say <- function(...) if (isTRUE(verbose)) message(...)

    ab <- mc_abundance_parts(abundance)
    md <- mc_metadata_parts(
        metadata, sample_col, subject_col, time_col, colnames(ab$mat)
    )

    if (identical(abundance_type, "raw")) {
        say("Transforming raw abundances to centred log-ratios.")
        ab$mat <- mc_clr(ab$mat, pseudocount)
    }

    design <- mc_build_design(ab, md, say)
    mc_report_design(design, say)
    design
}
