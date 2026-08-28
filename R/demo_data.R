# The bundled example, in the shape a study would supply it.
#
# taxa_demo is stored in the layout the fitting code works on, because that
# is what mc_prepare() and the benchmarking functions take. Presenting it
# that way to someone learning the package would teach the wrong habit, so
# this splits it into the two tables a study actually has: an abundance
# matrix with taxa in the row names, and a metadata table.

#' The bundled example as an abundance table and metadata
#'
#' @description
#' Returns `taxa_demo` in the form [mc_run()] expects: an abundance matrix
#' with taxa in the row names and samples in columns, and a metadata table
#' naming the subject and time behind each sample.
#'
#' Sample names are deliberately arbitrary, as a real sequencing run's would
#' be, to make clear that nothing is encoded in them.
#'
#' @return A list with `counts`, a numeric matrix, and `metadata`, a
#'   data.frame with columns `sample`, `subject` and `time`.
#'
#' @examples
#' demo <- mc_demo_data()
#'
#' demo$counts[1:3, 1:4]
#' head(demo$metadata)
#'
#' run <- suppressWarnings(mc_run(
#'     demo$counts, demo$metadata,
#'     sample_col = "sample", subject_col = "subject", time_col = "time",
#'     K = 1, verbose = FALSE
#' ))
#' run$missing
#'
#' @seealso [taxa_demo] for the stored form, used by [mc_prepare()].
#'
#' @export
mc_demo_data <- function() {
    # taxa_demo is lazy-loaded data rather than a namespace object, so it is
    # fetched into a private environment instead of read directly.
    env <- new.env(parent = emptyenv())
    utils::data("taxa_demo", package = "MicrobiomeCurves", envir = env)
    demo <- env$taxa_demo

    cols <- setdiff(names(demo), "OTU_ID")
    parsed <- mc_parse_cols(cols)

    mat <- as.matrix(demo[, cols, drop = FALSE])
    rownames(mat) <- demo$OTU_ID
    colnames(mat) <- sprintf("LIB_%03d", seq_along(cols))

    list(
        counts = mat,
        metadata = data.frame(
            sample = colnames(mat),
            subject = parsed$rep,
            time = parsed$time,
            stringsAsFactors = FALSE
        )
    )
}
