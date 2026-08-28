# Naming and describing samples the run had to create.
#
# A subject-timepoint with no sample has no name of the user's to carry, so
# one is built from the subject and the time. The name is made distinct from
# every real sample name, because the whole point is being able to tell an
# imputed sample from a measured one.

#' Name a column created for a subject-timepoint that had no sample
#'
#' @param cols Character vector of internal column names to name.
#' @param design The `mc_design` the run was built from.
#'
#' @return Character vector of names, one per entry of `cols`.
#'
#' @keywords internal
#' @noRd
mc_name_added <- function(cols, design) {
    parsed <- mc_parse_cols(cols)
    subject <- design$subjects[
        as.integer(sub("^s", "", parsed$rep))
    ]
    time <- design$times[parsed$time + 1L]

    base <- paste0(subject, "_", mc_time_label(design$axis, time))
    mc_make_distinct(base, design$map$sample)
}

#' Make generated names distinct from the real ones
#'
#' @description
#' A study could already contain a sample literally called `M01_7`, which
#' would make an imputed column indistinguishable from a measured one. Any
#' collision gets a suffix.
#'
#' @param proposed Character vector of generated names.
#' @param taken Character vector of names already in use.
#'
#' @return Character vector the same length as `proposed`.
#'
#' @keywords internal
#' @noRd
mc_make_distinct <- function(proposed, taken) {
    out <- proposed
    for (i in seq_along(out)) {
        candidate <- out[i]
        k <- 1L
        while (candidate %in% c(taken, out[-i])) {
            candidate <- paste0(proposed[i], "_imputed", if (k > 1L) k else "")
            k <- k + 1L
        }
        out[i] <- candidate
    }
    out
}

#' Extend the metadata with a row for each created sample
#'
#' @param design The `mc_design` the run was built from.
#' @param added Character vector of internal column names that were created.
#' @param labels Character vector of the final names of every output column.
#'
#' @return The design's metadata with one further row per created sample,
#'   carrying `imputed = TRUE`.
#'
#' @keywords internal
#' @noRd
mc_extend_metadata <- function(design, added, labels) {
    md <- design$metadata
    if (length(added) == 0) {
        return(md)
    }

    new_names <- mc_name_added(added, design)
    parsed <- mc_parse_cols(added)

    extra <- data.frame(
        sample = new_names,
        subject = design$subjects[as.integer(sub("^s", "", parsed$rep))],
        time = design$times[parsed$time + 1L],
        time_label = mc_time_label(
            design$axis, design$times[parsed$time + 1L]
        ),
        imputed = rep(TRUE, length(added)),
        stringsAsFactors = FALSE
    )

    out <- rbind(md, extra)
    rownames(out) <- NULL
    out
}
