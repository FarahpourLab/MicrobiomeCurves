#' @importFrom methods is new
#' @importFrom SummarizedExperiment assay assays assayNames assay<- colData
#' @importFrom SummarizedExperiment colData<- rowData
#' @importFrom S4Vectors metadata metadata<- DataFrame
NULL

#' Build safe subject.time column labels
#'
#' @description
#' Maps arbitrary subject identifiers and time values onto the
#' `<subject>.<time>` labels the table pipeline parses.
#'
#' @details
#' Subject identifiers in real data contain dots, spaces and other characters,
#' and time points are not always non-negative integers. Rather than constrain
#' what a `SummarizedExperiment` may hold, subjects are relabelled `s1`, `s2`,
#' ... and time points are replaced by their rank, `0`, `1`, ... The mapping is
#' returned so results can be translated back.
#'
#' @param subject Vector of subject identifiers, one per sample.
#' @param time Vector of time points, one per sample.
#'
#' @return A list with `label` (the synthetic column labels), `subjects` and
#'   `times` (the unique original values, in the order their codes were
#'   assigned).
#'
#' @keywords internal
#' @noRd
tti_encode_samples <- function(subject, time) {
    subjects <- unique(as.character(subject))
    times <- sort(unique(time))

    subj_code <- paste0("s", match(as.character(subject), subjects))
    time_code <- match(time, times) - 1L

    list(
        label = paste0(subj_code, ".", time_code),
        subjects = subjects,
        times = times
    )
}

#' Impute missing timepoints in a SummarizedExperiment
#'
#' @description
#' Method of [tti_run()] for `SummarizedExperiment` and
#' `TreeSummarizedExperiment` objects, the containers used across the
#' Bioconductor microbiome ecosystem.
#'
#' @details
#' The chosen assay is treated as a taxa-by-samples abundance matrix, expected
#' to hold CLR-transformed values. Subject and time point are read from
#' `colData`, so samples may be named however you like and time points need not
#' be consecutive integers.
#'
#' Missingness is detected exactly as for a plain table, in two forms: a sample
#' whose assay column is entirely `NA`, and a subject-timepoint combination
#' that has no column at all. Because a completely absent sample cannot be
#' represented in the returned object without a column to hold it, such samples
#' are appended, with `colData` filled in for `subject_col` and `time_col` and
#' `NA` elsewhere. **The returned object can therefore have more columns than
#' the input.** `metadata(x)$tti_run$added_samples` names them, and is
#' `character(0)` when none were added.
#'
#' The input assay is never modified. Imputed values are placed in a new assay,
#' so observed and imputed data stay distinguishable.
#'
#' @param dat A `SummarizedExperiment` or `TreeSummarizedExperiment`.
#' @param subject_col Character. Name of the `colData` column holding subject
#'   identifiers.
#' @param time_col Character. Name of the `colData` column holding time points.
#'   Values must be numeric or coercible to numeric; they need not be integers
#'   or evenly spaced.
#' @param assay_name Character or `NULL`. Which assay to impute. Defaults to
#'   the first assay.
#' @param name Character. Name of the assay to write the completed matrix into.
#'   Defaults to `"imputed"`.
#' @param ... Further arguments passed to the `data.frame` method of
#'   [tti_run()], such as `K`, `use_outliers`, `seed` and `verbose`.
#'
#' @return The input object with an additional assay named by `name`, holding
#'   the completed matrix, and with `metadata(x)$tti_run` set to a list
#'   carrying `missing`, `observed`, `n_failed` and `added_samples`. Columns
#'   are added if any subject-timepoint had no sample at all.
#'
#' @examples
#' # Build a small TreeSummarizedExperiment from the bundled demo table.
#' library(SummarizedExperiment)
#' library(TreeSummarizedExperiment)
#'
#' se <- tti_as_demo_se()
#' se
#'
#' out <- suppressWarnings(
#'     tti_run(se, subject_col = "subject", time_col = "timepoint",
#'         K = 1, verbose = FALSE)
#' )
#'
#' assayNames(out)
#' metadata(out)$tti_run$missing
#'
#' @seealso [tti_run()], [tti_as_demo_se()]
#'
#' @rdname tti_run
#' @export
setMethod("tti_run", "SummarizedExperiment", function(
    dat,
    subject_col = "subject",
    time_col = "timepoint",
    assay_name = NULL,
    name = "imputed",
    ...
) {
    cd <- SummarizedExperiment::colData(dat)
    spec <- tti_se_inputs(dat, cd, subject_col, time_col, assay_name)

    mat <- spec$mat
    enc <- spec$enc
    taxa <- spec$taxa

    tab <- data.frame(OTU_ID = taxa, stringsAsFactors = FALSE)
    tab[enc$label] <- as.data.frame(mat)

    run <- tti_run(tab, taxon_col = "OTU_ID", ...)

    completed <- run$completed
    out_labels <- setdiff(colnames(completed), "OTU_ID")
    added <- setdiff(out_labels, enc$label)

    filled <- as.matrix(completed[, c(enc$label, added), drop = FALSE])
    # taxa may be generated names when the object has no rownames, and the
    # assay setter rejects dimnames that differ from the receiving object's.
    # Mirroring the object keeps that case working.
    rownames(filled) <- rownames(dat)

    obj <- if (length(added) > 0) {
        tti_se_append_samples(dat, cd, enc, added, subject_col, time_col, taxa)
    } else {
        dat
    }

    colnames(filled) <- colnames(obj)
    SummarizedExperiment::assay(obj, name) <- filled

    S4Vectors::metadata(obj)$tti_run <- list(
        missing = run$missing,
        observed = run$observed,
        n_failed = run$n_failed,
        added_samples = added,
        assay_imputed = spec$assay_name
    )

    obj
})

#' A small TreeSummarizedExperiment for examples
#'
#' @description
#' Wraps [taxa_demo] in a `TreeSummarizedExperiment`, with `subject` and
#' `timepoint` columns in `colData`, so the container interface can be
#' demonstrated without constructing an object by hand.
#'
#' @details
#' The demo table is missing the sample `S04.2` entirely. That column is
#' created here holding `NA`, so the returned object is rectangular and every
#' subject-timepoint has a place. Two further samples, `S02.1` and `S07.4`,
#' are present but entirely `NA`.
#'
#' @return A `TreeSummarizedExperiment` with one assay, `"clr"`.
#'
#' @examples
#' se <- tti_as_demo_se()
#' dim(se)
#' head(SummarizedExperiment::colData(se))
#'
#' @seealso [taxa_demo], [tti_run()]
#'
#' @export
tti_as_demo_se <- function() {
    tab <- tti_demo_table()
    labels <- setdiff(colnames(tab), "OTU_ID")

    parts <- strsplit(labels, ".", fixed = TRUE)
    subject <- vapply(parts, `[`, character(1), 1)
    timepoint <- as.numeric(vapply(parts, `[`, character(1), 2))

    # Restore the fully absent sample as an all-NA column, so the object is
    # rectangular over subject x timepoint.
    grid <- expand.grid(
        subject = sort(unique(subject), method = "radix"),
        timepoint = sort(unique(timepoint)),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    grid$label <- paste0(grid$subject, ".", grid$timepoint)
    ord <- order(grid$subject, grid$timepoint, method = "radix")
    grid <- grid[ord, , drop = FALSE]

    mat <- matrix(
        NA_real_,
        nrow = nrow(tab), ncol = nrow(grid),
        dimnames = list(tab$OTU_ID, grid$label)
    )
    present <- intersect(grid$label, labels)
    mat[, present] <- as.matrix(tab[, present, drop = FALSE])

    TreeSummarizedExperiment::TreeSummarizedExperiment(
        assays = list(clr = mat),
        colData = S4Vectors::DataFrame(
            subject = grid$subject,
            timepoint = grid$timepoint,
            row.names = grid$label
        )
    )
}

#' Validate a container and pull out what the table pipeline needs
#'
#' @param dat A `SummarizedExperiment`.
#' @param cd Its `colData`.
#' @param subject_col Character name of the subject column.
#' @param time_col Character name of the time column.
#' @param assay_name Character name of the assay, or `NULL` for the
#'   first one.
#'
#' @return List with the assay matrix `mat`, the sample encoding `enc`,
#'   the taxon names `taxa`, and the resolved `assay_name`.
#'
#' @keywords internal
#' @noRd
tti_se_inputs <- function(dat, cd, subject_col, time_col, assay_name) {
    for (nm in c(subject_col, time_col)) {
        if (!(nm %in% colnames(cd))) {
            stop("colData has no column '", nm, "'")
        }
    }

    if (is.null(assay_name)) {
        assay_name <- SummarizedExperiment::assayNames(dat)[1]
        if (is.na(assay_name) || is.null(assay_name)) {
            stop("The object has no assay to impute.")
        }
    }
    if (!(assay_name %in% SummarizedExperiment::assayNames(dat))) {
        stop("No assay named '", assay_name, "' in the object.")
    }

    mat <- SummarizedExperiment::assay(dat, assay_name)
    if (!is.numeric(mat)) {
        stop("Assay '", assay_name, "' is not numeric.")
    }

    subject <- as.character(cd[[subject_col]])
    time <- tti_se_times(cd[[time_col]], time_col)

    enc <- tti_encode_samples(subject, time)
    if (anyDuplicated(enc$label) > 0) {
        stop(
            "Each subject may appear at most once per time point. ",
            "Duplicated subject-timepoint pairs were found."
        )
    }

    taxa <- rownames(mat)
    if (is.null(taxa)) {
        taxa <- paste0("feature", seq_len(nrow(mat)))
    }

    list(mat = mat, enc = enc, taxa = taxa, assay_name = assay_name)
}

#' Coerce a colData time column to numeric
#'
#' @param time_raw The column as stored.
#' @param time_col Character name of the column, used in errors.
#'
#' @return Numeric vector of time points.
#'
#' @keywords internal
#' @noRd
tti_se_times <- function(time_raw, time_col) {
    if (is.numeric(time_raw)) {
        time <- as.numeric(time_raw)
    } else {
        # Check before converting, so an invalid column gives a clear
        # error instead of a coercion warning.
        txt <- trimws(as.character(time_raw))
        num <- "^[+-]?(([0-9]+[.]?[0-9]*)|([.][0-9]+))([eE][+-]?[0-9]+)?$"
        if (!all(grepl(num, txt))) {
            stop(
                "colData column '", time_col,
                "' is not coercible to numeric."
            )
        }
        time <- as.numeric(txt)
    }
    if (anyNA(time)) {
        stop("colData column '", time_col,
             "' contains missing time points.")
    }
    time
}

#' Append columns for subject-timepoints that had no sample
#'
#' @param dat The input `SummarizedExperiment`.
#' @param cd Its `colData`.
#' @param enc Sample encoding from [tti_encode_samples()].
#' @param added Character vector of labels to append.
#' @param subject_col,time_col Column names in `colData`.
#' @param taxa Character vector of taxon names.
#'
#' @return A `SummarizedExperiment` with the extra columns.
#'
#' @keywords internal
#' @noRd
tti_se_append_samples <- function(dat, cd, enc, added, subject_col,
                                  time_col, taxa) {
    parts <- strsplit(added, ".", fixed = TRUE)
    add_subject <- enc$subjects[
        as.integer(sub("^s", "", vapply(parts, `[`, character(1), 1)))
    ]
    add_time <- enc$times[
        as.integer(vapply(parts, `[`, character(1), 2)) + 1L
    ]

    pad <- cd[rep(1L, length(added)), , drop = FALSE]
    rownames(pad) <- added
    for (nm in colnames(pad)) {
        pad[[nm]] <- rep(NA, length(added))
    }
    pad[[subject_col]] <- add_subject
    pad[[time_col]] <- add_time

    empty <- matrix(
        NA_real_,
        nrow = length(taxa), ncol = length(added),
        dimnames = list(taxa, added)
    )

    pad_assays <- lapply(
        SummarizedExperiment::assays(dat), function(a) cbind(a, empty)
    )

    obj <- SummarizedExperiment::SummarizedExperiment(
        assays = pad_assays,
        rowData = SummarizedExperiment::rowData(dat),
        colData = rbind(cd, pad)
    )
    S4Vectors::metadata(obj) <- S4Vectors::metadata(dat)
    obj
}
