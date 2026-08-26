# Shared fixtures.
#
# LazyData is off, as Bioconductor asks, so the demo table is loaded
# explicitly rather than being attached automatically.
utils::data("taxa_demo", package = "TaxaTimeImpute", envir = environment())

#
# fdapace emits "There is a time gap of at least ..." on sparse designs.
# This is expected for this kind of data, so the tests suppress it.

quiet_fit <- function(...) suppressWarnings(tti_fit(...))

quiet_run <- function(...) suppressWarnings(tti_run(..., verbose = FALSE))

# A single masked cell: enough to exercise the pipeline, fast enough for tests.
demo_prep <- function() {
    tti_prepare(
        dat = taxa_demo,
        taxon_col = "OTU_ID",
        mask_list = data.frame(rep = "S01", time = 3)
    )
}
