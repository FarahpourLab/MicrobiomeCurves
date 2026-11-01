# Shared fixtures.
#
# LazyData is off, as Bioconductor asks, so the demo table is loaded
# explicitly rather than being attached automatically.
utils::data("taxa_demo", package = "MicrobiomeCurves", envir = environment())

# The bundled table is the whole BONUS cohort: 833 taxa over 64 subjects,
# which takes hours to fit. The tests are about whether the pipeline does
# what it says, not about the cohort, so `taxa_demo` is cut down here to a
# corner of it and every test below works on the slice. The values, the
# column naming and the pattern of missing samples are the cohort's own, so
# the tests still run on real data, just less of it.
#
# The name is deliberately reused: tests express their expectations relative
# to `taxa_demo` -- nrow(taxa_demo) * 3, and so on -- so slicing it here
# keeps those expectations true without each test having to know the size.
demo_n_taxa <- 4L
demo_n_subjects <- 10L

local({
    cols <- setdiff(names(taxa_demo), "OTU_ID")
    subject <- sub("\\..*$", "", cols)

    subs <- sort(unique(subject))[seq_len(demo_n_subjects)]
    keep <- cols[subject %in% subs]

    # The taxa that vary most over this corner, so there is something to
    # fit. A taxon sitting at the detection floor in every sample is
    # constant, and a constant curve tests nothing.
    vals <- as.matrix(taxa_demo[, keep, drop = FALSE])
    v <- apply(vals, 1, stats::var)
    rows <- order(v, decreasing = TRUE)[seq_len(demo_n_taxa)]

    taxa_demo <<- taxa_demo[rows, c("OTU_ID", keep), drop = FALSE]
    rownames(taxa_demo) <<- NULL
})

# mc_demo_data() returns the whole cohort, which is hours of fitting. This
# is the same two-table shape over the slice above, with the cohort's own
# metadata rows for the samples the slice kept.
demo_small <- function() {
    # Built by cutting mc_demo_data()'s own result rather than by repeating
    # what it does, so the fixture cannot drift away from the function the
    # tests are about: the sample renaming, the column order and the
    # metadata join are all exercised as the user gets them.
    d <- mc_demo_data()
    keep_src <- setdiff(names(taxa_demo), "OTU_ID")
    j <- which(as.character(d$metadata$source_sample) %in% keep_src)

    meta <- d$metadata[j, , drop = FALSE]
    rownames(meta) <- NULL
    list(
        abundance = d$abundance[taxa_demo$OTU_ID, j, drop = FALSE],
        metadata = meta
    )
}

# mc_as_demo_se() reads the bundled table itself, so it does not see the
# slice above. The object is cut to the same corner here: the slice's taxa,
# and the columns belonging to the slice's subjects. It stays rectangular
# over subject x timepoint, which is the property the tests are about.
demo_se_small <- function() {
    se <- mc_as_demo_se()
    subs <- unique(sub("\\..*$", "", setdiff(names(taxa_demo), "OTU_ID")))
    keep <- SummarizedExperiment::colData(se)$subject %in% subs
    se[taxa_demo$OTU_ID, keep]
}

# A subject and time that the slice actually has, for mask_list fixtures.
demo_masked_cell <- local({
    cols <- setdiff(names(taxa_demo), "OTU_ID")
    first <- cols[1]
    data.frame(
        rep = sub("\\..*$", "", first),
        time = as.integer(sub("^.*\\.", "", first)),
        stringsAsFactors = FALSE
    )
})

#
# fdapace emits "There is a time gap of at least ..." on sparse designs.
# This is expected for this kind of data, so the tests suppress it.

quiet_fit <- function(...) suppressWarnings(mc_fit(...))

# mc_run() takes an abundance table plus metadata, so a fixture stored in
# the internal layout is converted first. Numbers are unaffected: the
# conversion only renames columns.
quiet_run <- function(...) {
    suppressWarnings(run_wide_as_meta(..., verbose = FALSE))
}

# A single masked cell: enough to exercise the pipeline, fast enough for tests.
demo_prep <- function() {
    mc_prepare(
        dat = taxa_demo,
        taxon_col = "OTU_ID",
        mask_list = demo_masked_cell
    )
}
