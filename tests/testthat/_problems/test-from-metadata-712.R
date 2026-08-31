# Extracted from test-from-metadata.R:712

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "MicrobiomeCurves", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
study <- function(drop = NULL, blank = NULL) {
    subjects <- paste0("M", sprintf("%02d", 1:5))
    days <- c(0, 7, 14)

    meta <- expand.grid(
        subject = subjects, day = days,
        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
    meta <- meta[order(meta$subject, meta$day), ]
    meta$sample <- sprintf("RUN_%04d", seq_len(nrow(meta)))
    meta <- meta[, c("sample", "subject", "day")]
    rownames(meta) <- NULL

    pick <- function(sel) {
        meta$sample[meta$subject == sel[[1]] & meta$day == sel[[2]]]
    }
    dropped <- if (!is.null(drop)) pick(drop) else NULL
    blanked <- if (!is.null(blank)) pick(blank) else NULL

    if (!is.null(dropped)) meta <- meta[meta$sample != dropped, ]
    rownames(meta) <- NULL

    taxa <- paste0("Taxon", sprintf("%02d", 1:4))
    set.seed(42)
    m <- matrix(
        rnorm(length(taxa) * nrow(meta)),
        nrow = length(taxa), dimnames = list(taxa, meta$sample)
    )
    if (!is.null(blanked)) m[, blanked] <- NA_real_

    list(
        abundance = m, metadata = meta, taxa = taxa,
        subjects = subjects, days = days,
        dropped = dropped, blanked = blanked
    )
}
build <- function(s, ...) {
    mc_from_metadata(
        s$abundance, s$metadata,
        sample_col = "sample", subject_col = "subject", time_col = "day",
        verbose = FALSE, ...
    )
}

# test -------------------------------------------------------------------------
set.seed(4)
subs <- paste0("SUB", sprintf("%02d", 1:6))
days <- c(0, 7, 14, 21)
meta <- expand.grid(
        SubjectID = subs, Day = days,
        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
meta <- meta[order(meta$SubjectID, meta$Day), ]
meta$SampleID <- sprintf("S%03d", seq_len(nrow(meta)))
meta <- meta[, c("SampleID", "SubjectID", "Day")]
meta <- meta[!(meta$SubjectID == "SUB06" & meta$Day > 0), ]
taxa <- paste0("Genus_", LETTERS[1:3])
counts <- matrix(
        rnorm(length(taxa) * nrow(meta), mean = 2),
        nrow = length(taxa), dimnames = list(taxa, meta$SampleID)
    )
warned <- character(0)
run <- withCallingHandlers(
        mc_run(
            counts, meta,
            sample_col = "SampleID", subject_col = "SubjectID",
            time_col = "Day", K = 1, verbose = FALSE
        ),
        warning = function(w) {
            warned <<- c(warned, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
