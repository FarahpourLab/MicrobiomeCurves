# ============================================================
# Generates the bundled example dataset `taxa_demo`.
#
# Run from the package root:
#     Rscript data-raw/make_taxa_demo.R
#
# Writes:
#     data/taxa_demo.rda          lazy-loaded object
#     inst/extdata/taxa_demo.csv  same table as a raw CSV
#
# The data are SIMULATED, not from any real cohort. The design
# mirrors the BONUS cystic-fibrosis cohort used in the paper:
# a handful of taxa, a modest number of subjects, 7 time points,
# and a few whole samples missing.
# ============================================================

set.seed(123)

n_taxa <- 12
subjects <- sprintf("S%02d", 1:10)
times <- 0:6

# Each taxon follows a smooth curve over time; each subject is
# offset by a random shift. Values are on a CLR-like scale, so
# they are centred around zero and may be negative.
amplitude <- runif(n_taxa, 0.5, 2.0)
phase <- runif(n_taxa, 0, pi)
shift <- rnorm(length(subjects), 0, 0.5)

taxa_demo <- data.frame(
    OTU_ID = sprintf("Taxon%02d", seq_len(n_taxa)),
    stringsAsFactors = FALSE
)

for (i in seq_along(subjects)) {
    for (tt in times) {
        taxa_demo[[paste0(subjects[i], ".", tt)]] <-
            amplitude * sin(tt / 2 + phase) +
            shift[i] +
            rnorm(n_taxa, 0, 0.3)
    }
}

# ---- introduce missing whole samples --------------------------------
# Two forms are represented, because tti_detect_missing() recognises both:
#
#   S04.2  -> column removed entirely  ("absent_column")
#   S07.4  -> column present but all NA ("all_na")
#   S02.1  -> column present but all NA ("all_na")

taxa_demo[["S04.2"]] <- NULL
taxa_demo[["S07.4"]] <- NA_real_
taxa_demo[["S02.1"]] <- NA_real_

stopifnot(
    nrow(taxa_demo) == n_taxa,
    !("S04.2" %in% names(taxa_demo)),
    all(is.na(taxa_demo[["S07.4"]])),
    all(is.na(taxa_demo[["S02.1"]]))
)

dir.create("data", showWarnings = FALSE)
dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)

save(taxa_demo, file = "data/taxa_demo.rda", compress = "bzip2", version = 2)
write.csv(taxa_demo, "inst/extdata/taxa_demo.csv", row.names = FALSE)

cat("taxa_demo:", nrow(taxa_demo), "taxa x",
    ncol(taxa_demo) - 1, "sample columns\n")
cat("wrote data/taxa_demo.rda and inst/extdata/taxa_demo.csv\n")
