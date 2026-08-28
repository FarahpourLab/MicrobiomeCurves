# Render the workflow vignette to HTML and PDF in the package root.
#
# These are reading copies, for sharing outside the package. They are not
# vignettes: only the HTML vignette ships, because a PDF vignette would make
# every build machine need a LaTeX engine, and none of the three used by CI
# has one.
#
# Run from the package root:
#     Rscript tools/render-docs.R

src <- "vignettes/MicrobiomeCurves-workflow.Rmd"
if (!file.exists(src)) {
    stop("Run this from the package root; ", src, " was not found.")
}

for (pkg in c("rmarkdown", "BiocStyle", "callr")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("Package '", pkg, "' is needed to render the documents.")
    }
}

# LaTeX writes its intermediates beside the source, and R CMD check calls
# anything left in vignettes/ a leftover.
tidy_up <- function() {
    unlink(list.files(
        "vignettes",
        pattern = "[.](log|tex|aux|out|toc|knit[.]md|synctex[.]gz)$",
        full.names = TRUE
    ))
}
on.exit(tidy_up(), add = TRUE)
tidy_up()

# Each format renders in its own process. Both write the same intermediate
# file beside the source, so sharing a session makes the second one fail.
render_one <- function(src, fmt) {
    callr::r(
        function(src, fmt) {
            out <- if (fmt == "html") {
                BiocStyle::html_document(toc_float = TRUE)
            } else {
                BiocStyle::pdf_document(toc = TRUE)
            }
            rmarkdown::render(src,
                output_format = out, output_dir = ".", quiet = TRUE
            )
        },
        args = list(src = src, fmt = fmt),
        show = FALSE
    )
}

made <- character(0)
for (fmt in c("html", "pdf")) {
    message("Rendering ", toupper(fmt), " ...")
    made <- c(made, render_one(src, fmt))
    tidy_up()
}

for (f in made) {
    message(sprintf("  %-34s %8.0f bytes", basename(f), file.info(f)$size))
}
