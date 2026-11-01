# Builds the bundled example data from the BONUS cystic-fibrosis infant
# cohort, so that installing the package is enough to have it.
#
#     Rscript data-raw/make_taxa_demo.R
#
# Writes
#     data/taxa_demo.rda            the abundance table, lazy-loaded
#     data/taxa_demo_meta.rda       the sample metadata
#     inst/extdata/taxa_demo.csv    the same abundance table as a raw CSV
#
# The two inputs are produced outside the package by
# 03.processed/05.package_input/make_package_input.R, which applies one
# missingness mask and drops the cluster column. They are not part of the
# package sources, so this script only runs where the cohort is available.
# It is kept here so the bundled data can be traced back to the files it
# came from.
#
# Two things are done to the cohort, both of which the documentation states:
#
#   The columns are named "<subject>.<time>", which the cohort already uses
#   and which mc_detect_missing() parses. A subject-timepoint that the mask
#   removed simply has no column, which is how a sample that was never
#   collected reaches the package.
#
#   Subjects with fewer than two observed time points are dropped. mc_run()
#   refuses them outright, since a single observation carries almost nothing
#   of its own, and example data that cannot be run is of no use. The MAR
#   0.3 masks each leave one such subject.

in_dir <- file.path(
    "F:/Dropbox/Personal/UDE/Projects/CANDi/01.benchmarking_datasets",
    "bonus_cf_infants/03.processed/05.package_input/MAR/0.3/rep_1"
)
ab_in <- file.path(in_dir, "clr_species_abundance.csv")
md_in <- file.path(in_dir, "sample_metadata.csv")

for (f in c(ab_in, md_in)) {
    if (!file.exists(f)) {
        stop(
            "there is no file at\n  ", f,
            "\nRun 05.package_input/make_package_input.R first.",
            call. = FALSE
        )
    }
}

ab <- read.csv(ab_in, check.names = FALSE, stringsAsFactors = FALSE)
md <- read.csv(md_in, check.names = FALSE, stringsAsFactors = FALSE)
cat("read:", nrow(ab), "taxa x", ncol(ab) - 1, "samples\n")

# ---- drop subjects with a single observed time point --------------------
per <- table(md$subject)
thin <- names(per)[per < 2]
if (length(thin) > 0) {
    cat("dropping", length(thin), "subject(s) with one observation:",
        paste(thin, collapse = ", "), "\n")
    md <- md[!(md$subject %in% thin), , drop = FALSE]
    ab <- ab[, c("taxon", md$sample), drop = FALSE]
}

# ---- the stored form: taxon identifiers in OTU_ID -----------------------
names(ab)[names(ab) == "taxon"] <- "OTU_ID"
taxa_demo <- ab
rownames(taxa_demo) <- NULL

taxa_demo_meta <- md
rownames(taxa_demo_meta) <- NULL

samples <- setdiff(names(taxa_demo), "OTU_ID")
n_sub <- length(unique(taxa_demo_meta$subject))
n_time <- length(unique(taxa_demo_meta$time))

stopifnot(
    identical(samples, as.character(taxa_demo_meta$sample)),
    !anyDuplicated(taxa_demo$OTU_ID),
    all(vapply(taxa_demo[samples], is.numeric, logical(1))),
    all(table(taxa_demo_meta$subject) >= 2)
)

cat("\ntaxa_demo      :", nrow(taxa_demo), "taxa x", length(samples),
    "sample columns\n")
cat("taxa_demo_meta :", nrow(taxa_demo_meta), "rows x",
    ncol(taxa_demo_meta), "columns\n")
cat("subjects       :", n_sub, "\n")
cat("time points    :", n_time, "\n")
cat("grid           :", n_sub * n_time, "cells, of which",
    n_sub * n_time - length(samples), "have no column\n")
cat("NA in values   :", sum(is.na(taxa_demo[samples])), "\n")

dir.create("data", showWarnings = FALSE)
dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)

save(taxa_demo,
    file = "data/taxa_demo.rda", compress = "xz", version = 2
)
save(taxa_demo_meta,
    file = "data/taxa_demo_meta.rda", compress = "xz", version = 2
)
# Gzipped, because the table is 4.9 MB as plain text and a package has a
# size limit to respect. read.csv() and read.table() open a .gz directly,
# so nothing about the example changes.
unlink("inst/extdata/taxa_demo.csv")
con <- gzfile("inst/extdata/taxa_demo.csv.gz", "wt")
utils::write.csv(taxa_demo, con, row.names = FALSE)
close(con)

con <- gzfile("inst/extdata/taxa_demo_meta.csv.gz", "wt")
utils::write.csv(taxa_demo_meta, con, row.names = FALSE)
close(con)

cat("\nwrote\n")
for (f in c("data/taxa_demo.rda", "data/taxa_demo_meta.rda",
    "inst/extdata/taxa_demo.csv.gz", "inst/extdata/taxa_demo_meta.csv.gz")) {
    cat(sprintf("  %-34s %7.2f MB\n", f, file.size(f) / 1024^2))
}
