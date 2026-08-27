# Validation and assembly for tti_from_metadata().
#
# Each check reports its own cause. A mismatch between the abundance columns
# and the metadata rows is the commonest way this goes wrong, so it is named
# explicitly in both directions rather than surfacing later as an empty
# design.

#' Split an abundance table into a numeric matrix and taxon names
#'
#' @param abundance A data.frame or matrix, taxa in rows, samples in columns.
#' @param taxon_col Character name of the taxon column, or `NULL` to use row
#'   names.
#'
#' @return List with `mat`, the numeric matrix, and `taxa`.
#'
#' @keywords internal
#' @noRd
tti_abundance_parts <- function(abundance, taxon_col) {
    if (!is.data.frame(abundance) && !is.matrix(abundance)) {
        stop(
            "abundance must be a data.frame or a matrix with taxa in rows ",
            "and samples in columns, not ", class(abundance)[1], ".",
            call. = FALSE
        )
    }

    df <- as.data.frame(abundance, stringsAsFactors = FALSE)
    taxa <- tti_abundance_taxa(df, taxon_col, abundance)
    if (!is.null(taxon_col)) {
        df <- df[, setdiff(names(df), taxon_col), drop = FALSE]
    }

    if (ncol(df) == 0) {
        stop("abundance has no sample columns.", call. = FALSE)
    }
    bad <- names(df)[!vapply(df, is.numeric, logical(1))]
    if (length(bad) > 0) {
        stop(
            "These abundance columns are not numeric: ", tti_fmt_some(bad),
            ". Every column other than the taxon column is read as a ",
            "sample. If one of these is an annotation, remove it before ",
            "calling, or name it with taxon_col.",
            call. = FALSE
        )
    }

    list(mat = as.matrix(df), taxa = taxa)
}

#' Resolve taxon identifiers for an abundance table
#'
#' @param df The abundance table as a data.frame.
#' @param taxon_col Character name of the taxon column, or `NULL`.
#' @param original The object as supplied, used for its row names.
#'
#' @return Character vector of taxon identifiers.
#'
#' @keywords internal
#' @noRd
tti_abundance_taxa <- function(df, taxon_col, original) {
    if (!is.null(taxon_col)) {
        if (!(taxon_col %in% names(df))) {
            stop(
                "taxon_col '", taxon_col, "' is not a column of abundance. ",
                "Columns are: ", tti_fmt_some(names(df)),
                call. = FALSE
            )
        }
        return(as.character(df[[taxon_col]]))
    }

    rn <- rownames(original)
    if (is.null(rn) || all(rn == seq_len(nrow(df)))) {
        stop(
            "abundance has no row names to use as taxon identifiers. ",
            "Either set row names, or name the column holding them with ",
            "taxon_col.",
            call. = FALSE
        )
    }
    as.character(rn)
}

#' Validate the metadata and align it to the abundance columns
#'
#' @param metadata The user's metadata data.frame.
#' @param sample_col,subject_col,time_col Column names within it.
#' @param ab_samples Character vector of abundance column names.
#'
#' @return List with the aligned `sample`, `subject` and `time` vectors, and
#'   `absent`, the samples described in metadata but missing from the
#'   abundance table.
#'
#' @keywords internal
#' @noRd
tti_metadata_parts <- function(metadata, sample_col, subject_col, time_col,
                               ab_samples) {
    if (!is.data.frame(metadata)) {
        stop(
            "metadata must be a data.frame with one row per sample, not ",
            class(metadata)[1], ".",
            call. = FALSE
        )
    }

    tti_check_meta_cols(metadata, sample_col, subject_col, time_col)

    sample <- as.character(metadata[[sample_col]])
    subject <- as.character(metadata[[subject_col]])
    time <- tti_se_times(metadata[[time_col]], time_col)

    tti_check_meta_values(sample, subject, sample_col, subject_col)
    tti_check_sample_overlap(sample, ab_samples, sample_col)

    list(
        sample = sample, subject = subject, time = time,
        absent = setdiff(sample, ab_samples)
    )
}

#' Check the three named metadata columns exist and are distinct
#'
#' @param metadata The user's metadata data.frame.
#' @param sample_col,subject_col,time_col Column names within it.
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
#' @noRd
tti_check_meta_cols <- function(metadata, sample_col, subject_col, time_col) {
    named <- c(
        sample_col = sample_col, subject_col = subject_col,
        time_col = time_col
    )
    missing <- named[!(named %in% names(metadata))]
    if (length(missing) > 0) {
        quoted <- tti_fmt_some(paste0("'", missing, "'"))
        by <- tti_fmt_some(names(missing))
        have <- tti_fmt_some(names(metadata))

        stop(
            "metadata has no column ", quoted,
            " (named by ", by, "). Its columns are: ", have, ".",
            call. = FALSE
        )
    }
    if (anyDuplicated(named) > 0) {
        stop(
            "sample_col, subject_col and time_col must name three different ",
            "columns; got ", tti_fmt_some(named), ".",
            call. = FALSE
        )
    }
    invisible(NULL)
}

#' Check the metadata identifies each sample once
#'
#' @param sample,subject Character vectors from the metadata.
#' @param sample_col,subject_col Their column names, used in errors.
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
#' @noRd
tti_check_meta_values <- function(sample, subject, sample_col, subject_col) {
    if (anyNA(sample) || anyNA(subject)) {
        stop(
            "metadata columns '", sample_col, "' and '", subject_col,
            "' must not contain NA.",
            call. = FALSE
        )
    }

    dup <- unique(sample[duplicated(sample)])
    if (length(dup) > 0) {
        stop(
            "metadata lists these samples more than once: ",
            tti_fmt_some(dup), ". Each sample needs exactly one row.",
            call. = FALSE
        )
    }
    invisible(NULL)
}

#' Compare the metadata's samples with the abundance table's columns
#'
#' @param meta_samples Character vector of sample names in the metadata.
#' @param ab_samples Character vector of abundance column names.
#' @param sample_col Name of the metadata column, used in the error.
#'
#' @return `NULL`, invisibly. Errors when the two share nothing, or when the
#'   abundance table holds samples the metadata does not describe.
#'
#' @keywords internal
#' @noRd
tti_check_sample_overlap <- function(meta_samples, ab_samples, sample_col) {
    shared <- intersect(meta_samples, ab_samples)
    if (length(shared) == 0) {
        stop(
            "No sample name is shared between metadata['", sample_col,
            "'] and the columns of abundance. Metadata has ",
            tti_fmt_some(meta_samples), "; abundance has ",
            tti_fmt_some(ab_samples),
            ". Check that sample_col names the right column.",
            call. = FALSE
        )
    }

    undescribed <- setdiff(ab_samples, meta_samples)
    if (length(undescribed) > 0) {
        stop(
            length(undescribed), " abundance column(s) are not described in ",
            "the metadata: ", tti_fmt_some(undescribed),
            ". Every sample needs a row saying which subject and time it ",
            "belongs to. Add them, or drop those columns.",
            call. = FALSE
        )
    }
    invisible(NULL)
}
