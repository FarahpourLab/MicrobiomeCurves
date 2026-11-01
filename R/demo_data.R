# The bundled example, in the shape a study would supply it.
#
# taxa_demo is stored in the layout the fitting code works on, because that
# is what mc_prepare() and the benchmarking functions take: one table whose
# column names carry the subject and the time. Presenting it that way to
# someone learning the package would teach the wrong habit, so this splits
# it into the two tables a study actually has, an abundance matrix and a
# metadata table, and returns the cohort's own metadata rather than one
# reconstructed from the column names.

#' The bundled example as an abundance table and metadata
#'
#' @description
#' Returns the bundled BONUS data in the form [mc_run()] expects: an
#' abundance matrix with taxa in the row names and samples in columns, and
#' the cohort's metadata, one row per sample.
#'
#' @details
#' The values are centred log ratios, so [mc_run()] takes them with
#' `abundance_type = "clr"` and [mc_trends()] with its default
#' `normalization = "NONE"`.
#'
#' The metadata carries `sample`, `subject` and `time`, which every entry
#' point needs, and ten columns of the cohort's own measurements that can be
#' used as `covariates` or as `group_col`: `month`, `gender`, `panc_insuf`,
#' `genotype`, `antibiotics`, `breastmilk`, `formula`, `table_food`,
#' `weight` and `length`.
#'
#' `time` is the visit index, 0 to 6, which is the order the model works in.
#' `month` is the infant's age in months at that visit, and is a covariate
#' rather than the time axis.
#'
#' The whole table is 833 taxa over 64 subjects, and fitting all of it takes
#' hours. `n_taxa` and `n_subjects` return a corner of it instead, which is
#' the same data and enough to see how the functions behave.
#'
#' @param n_taxa Number of taxa to return, or `NULL` for all 833. The taxa
#'   that vary most over the samples returned are taken first, because most
#'   of this cohort sits at the detection floor in most samples and a
#'   constant curve shows nothing.
#' @param n_subjects Number of subjects to return, or `NULL` for all 64.
#'   Taken in order, so the same number always gives the same subjects.
#'
#' @return A list with `abundance`, a numeric matrix of taxa by samples, and
#'   `metadata`, a data.frame with one row per sample and 14 columns. With
#'   both arguments `NULL` that is 833 taxa and 319 samples.
#'
#' @examples
#' demo <- mc_demo_data()
#'
#' dim(demo$abundance)
#' demo$abundance[1:2, 1:4]
#' head(demo$metadata[, c("sample", "subject", "time", "month", "gender")])
#'
#' # The design is every subject crossed with every time point. The samples
#' # that do not exist are the difference, and are what mc_run() fills.
#' n_sub <- length(unique(demo$metadata$subject))
#' n_time <- length(unique(demo$metadata$time))
#' c(
#'     design = n_sub * n_time,
#'     collected = nrow(demo$metadata),
#'     missing = n_sub * n_time - nrow(demo$metadata)
#' )
#'
#' # a corner of it, for a quick look
#' small <- mc_demo_data(n_taxa = 3, n_subjects = 6)
#' dim(small$abundance)
#'
#' @seealso [taxa_demo] for the stored form, used by [mc_prepare()], and
#'   [taxa_demo_meta] for the metadata on its own.
#'
#' @export
mc_demo_data <- function(n_taxa = NULL, n_subjects = NULL) {
    # taxa_demo is lazy-loaded data rather than a namespace object, so it is
    # fetched into a private environment instead of read directly.
    env <- new.env(parent = emptyenv())
    utils::data(
        "taxa_demo", "taxa_demo_meta",
        package = "MicrobiomeCurves", envir = env
    )
    demo <- env$taxa_demo
    meta <- env$taxa_demo_meta

    cols <- setdiff(names(demo), "OTU_ID")
    mat <- as.matrix(demo[, cols, drop = FALSE])
    rownames(mat) <- demo$OTU_ID

    meta <- meta[match(cols, as.character(meta$sample)), , drop = FALSE]

    # Subjects first, then taxa, because which taxa vary most depends on
    # which samples are being looked at.
    if (!is.null(n_subjects)) {
        keep <- unique(as.character(meta$subject))
        keep <- keep[seq_len(min(length(keep), max(1L, as.integer(n_subjects))))]
        j <- which(as.character(meta$subject) %in% keep)
        meta <- meta[j, , drop = FALSE]
        mat <- mat[, j, drop = FALSE]
    }
    if (!is.null(n_taxa)) {
        v <- apply(mat, 1, stats::var, na.rm = TRUE)
        v[is.na(v)] <- -Inf
        n <- min(nrow(mat), max(1L, as.integer(n_taxa)))
        mat <- mat[order(v, decreasing = TRUE)[seq_len(n)], , drop = FALSE]
    }

    cols <- colnames(mat)

    # The stored table encodes the subject and the time in each column name,
    # because that is what the benchmarking route parses. This route does
    # not parse anything: it reads the subject and the time from the
    # metadata. So the samples are renamed to something arbitrary, as a real
    # sequencing run's identifiers would be, to make clear that nothing here
    # depends on the shape of a sample name. The cohort's own identifier is
    # kept alongside, so a sample can still be traced back to taxa_demo.
    meta$source_sample <- as.character(meta$sample)
    colnames(mat) <- sprintf("LIB_%03d", seq_along(cols))
    meta$sample <- colnames(mat)
    rownames(meta) <- NULL

    front <- c("sample", "subject", "time")
    meta <- meta[, c(front, setdiff(names(meta), front)), drop = FALSE]

    list(abundance = mat, metadata = meta)
}
