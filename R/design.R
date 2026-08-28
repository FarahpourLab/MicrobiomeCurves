# Assembling a mc_design, and reporting what it is missing.
#
# Two kinds of gap are distinguished, because they mean different things to
# whoever prepared the data:
#
#   no_data        the metadata lists the sample and the abundance table has
#                  a column for it, but every taxon in it is NA. The sample
#                  was planned and its column exists, but carries nothing.
#
#   absent_sample  no sample exists for that subject at that time. This is
#                  found from the grid of subjects against observed time
#                  points, not from the metadata, which only lists what was
#                  collected.
#
# A column imputed for the second kind has no name of the user's to carry, so
# one is built from the subject and the time and recorded in the returned
# metadata, flagged, so imputed samples stay distinguishable from real ones.

#' Assemble the internal table and the mapping back to sample names
#'
#' @param ab List from [mc_abundance_parts()].
#' @param md List from [mc_metadata_parts()].
#' @param say Function used to report progress.
#'
#' @return An object of class `mc_design`.
#'
#' @keywords internal
#' @noRd
mc_build_design <- function(ab, md, say) {
    ord <- match(colnames(ab$mat), md$sample)
    subject <- md$subject[ord]
    time <- md$time[ord]

    mc_check_design(unique(subject))
    mc_check_grid(subject, time)

    enc <- mc_encode_samples(subject, time)

    tab <- data.frame(OTU_ID = ab$taxa, stringsAsFactors = FALSE)
    tab[enc$label] <- as.data.frame(ab$mat)

    map <- data.frame(
        column = enc$label,
        sample = colnames(ab$mat),
        subject = subject,
        time = time,
        time_label = mc_time_label(md$axis, time),
        imputed = FALSE,
        stringsAsFactors = FALSE
    )

    # A column of all NA is present but carries nothing, so it counts as a
    # gap rather than an observation in both places below.
    empty <- colSums(!is.na(ab$mat)) == 0

    structure(
        list(
            table = tab,
            map = map,
            metadata = map[, c(
                "sample", "subject", "time", "time_label", "imputed"
            )],
            missing = mc_design_missing(empty, map, enc, md$axis),
            observed = mc_design_observed(map, empty),
            subjects = enc$subjects,
            times = enc$times,
            axis = md$axis,
            taxon_col = "OTU_ID"
        ),
        class = "mc_design"
    )
}

#' Refuse a design that measures a subject twice at one time
#'
#' @param subject,time Aligned vectors, one entry per sample.
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
#' @noRd
mc_check_grid <- function(subject, time) {
    key <- paste(subject, time, sep = "\r")
    dup <- unique(key[duplicated(key)])
    if (length(dup) == 0) {
        return(invisible(NULL))
    }

    parts <- do.call(rbind, strsplit(dup, "\r", fixed = TRUE))
    listed <- mc_fmt_some(paste0(parts[, 1], " at ", parts[, 2]))

    stop(
        "These subject-timepoint pairs appear more than once: ", listed,
        ". Each subject may be sampled at most once per time point; ",
        "technical replicates must be combined first.",
        call. = FALSE
    )
}

#' Find every subject-timepoint the design does not cover
#'
#' @param empty Logical vector, `TRUE` where a sample column holds no data.
#' @param map The column map.
#' @param enc List from [mc_encode_samples()].
#' @param axis The time axis, so the label the user wrote is carried too.
#'
#' @return Data frame with `subject`, `time`, `time_label`, `reason` and
#'   `sample`, ordered by subject then time.
#'
#' @keywords internal
#' @noRd
mc_design_missing <- function(empty, map, enc, axis) {
    no_data <- data.frame(
        subject = map$subject[empty], time = map$time[empty],
        reason = rep("no_data", sum(empty)),
        sample = map$sample[empty],
        stringsAsFactors = FALSE
    )

    grid <- expand.grid(
        subject = enc$subjects, time = enc$times,
        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
    have <- paste(map$subject, map$time, sep = "\r")
    gap <- !(paste(grid$subject, grid$time, sep = "\r") %in% have)

    absent <- data.frame(
        subject = grid$subject[gap], time = grid$time[gap],
        reason = rep("absent_sample", sum(gap)),
        sample = rep(NA_character_, sum(gap)),
        stringsAsFactors = FALSE
    )

    out <- rbind(no_data, absent)
    out$time_label <- mc_time_label(axis, out$time)
    out <- out[, c("subject", "time", "time_label", "reason", "sample")]
    if (nrow(out) == 0) {
        return(out)
    }
    out <- out[order(out$subject, out$time, method = "radix"), , drop = FALSE]
    rownames(out) <- NULL
    out
}

#' Count the time points each subject actually carries data for
#'
#' @param map The column map.
#' @param empty Logical vector, `TRUE` where a sample column holds no data.
#'
#' @return Data frame with `subject` and `n_observed`, one row per subject
#'   including any whose samples are all empty.
#'
#' @keywords internal
#' @noRd
mc_design_observed <- function(map, empty) {
    subjects <- unique(map$subject)
    n <- vapply(
        subjects,
        function(s) sum(map$subject == s & !empty),
        numeric(1)
    )
    data.frame(
        subject = subjects,
        n_observed = as.numeric(n),
        stringsAsFactors = FALSE
    )
}
