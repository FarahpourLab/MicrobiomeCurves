#' Longitudinal microbiome table from the BONUS infant cohort
#'
#' @description
#' The bundled example, in the layout the package's benchmarking route
#' works on: one table whose column names carry the subject and the time,
#' with whole samples missing. Installing the package is enough to have it,
#' so the workflow can be run without any external file.
#'
#' @details
#' The data are from the BONUS cystic-fibrosis infant cohort, at species
#' level, as centred log ratios. They are real measurements, not simulated.
#'
#' Columns are named `"<subject>.<time>"`, which is what
#' \code{\link{mc_detect_missing}} parses. A subject-timepoint with no
#' column is a sample that does not exist, which is the missingness the
#' package imputes. Of the 448 subject-timepoints in the design, 319 have a
#' column and 129 do not.
#'
#' Two things were done to the cohort before bundling, both deliberate:
#'
#' \itemize{
#'   \item One missingness mask was applied, so the table has gaps to
#'     impute. The mask itself is not part of the data; a masked sample
#'     simply has no column.
#'   \item Subjects with fewer than two observed time points were dropped.
#'     \code{\link{mc_run}} refuses such a subject, since one observation
#'     carries almost nothing of its own, and example data that cannot be
#'     run would be of no use.
#' }
#'
#' Every value is present: no column is partly or wholly \code{NA}. The
#' missingness is entirely whole samples that were never collected.
#'
#' The same table is also shipped as a raw CSV, gzipped, so the more
#' realistic \code{read.csv()} entry point can be demonstrated:
#'
#' \preformatted{
#' path <- system.file("extdata", "taxa_demo.csv.gz", package = "MicrobiomeCurves")
#' dat  <- read.csv(path, check.names = FALSE)
#' }
#'
#' @note
#' Fitting this dataset emits repeated \code{"There is a time gap of at
#' least ..."} warnings. These come from \code{fdapace::FPCA}, not from this
#' package, and are expected for sparse longitudinal designs. They are
#' notices about the sampling design, not errors, and the imputation is
#' unaffected.
#'
#' Fitting all 833 taxa takes hours. Subset the table, or pass \code{taxa}
#' to \code{\link{mc_trends}}, for a quick look.
#'
#' @format A data frame with 833 rows (taxa) and 320 columns.
#'
#' The first column, \code{OTU_ID}, is a character taxon identifier: a
#' full lineage from kingdom to species, such as
#' \code{"k__Bacteria|p__Firmicutes|...|s__Ezakiella_SGB6726"}.
#'
#' The remaining 319 columns are numeric centred log ratios, one per
#' subject-timepoint, named \code{"<subject>.<time>"} - 64 subjects crossed
#' with time points \code{0} to \code{6}, of which 319 of the 448
#' combinations are present.
#'
#' @source
#' The BONUS cystic-fibrosis infant cohort. Prepared by
#' \code{data-raw/make_taxa_demo.R} in the package sources, from the
#' processed cohort tables and one missingness mask.
#'
#' @usage data(taxa_demo)
#'
#' @seealso \code{\link{mc_demo_data}} for the same data as an abundance
#'   matrix and metadata, \code{\link{taxa_demo_meta}} for the metadata,
#'   \code{\link{mc_run}}, \code{\link{mc_detect_missing}}
#'
#' @examples
#' data(taxa_demo)
#'
#' taxa_demo[1:3, 1:4]
#'
#' # what is missing, and why
#' info <- mc_detect_missing(taxa_demo, taxon_col = "OTU_ID")
#' head(info$missing)
#'
#' # observations remaining per subject
#' head(info$observed)
"taxa_demo"


#' Sample metadata for the bundled BONUS table
#'
#' @description
#' One row per sample of \code{\link{taxa_demo}}, naming the subject and the
#' time behind each, with the cohort's own measurements alongside.
#'
#' @details
#' \code{time} is the visit index, 0 to 6, which is the order the model
#' works in. \code{month} is the infant's age in months at that visit, and
#' is a covariate rather than the time axis: the two are not the same, since
#' visits were not evenly spaced for every infant.
#'
#' The remaining columns can be passed as \code{covariates} to
#' \code{\link{mc_trends}} or as \code{group_col} to
#' \code{\link{mc_compare}}.
#'
#' @format A data frame with 319 rows and 13 columns:
#' \describe{
#'   \item{sample}{Sample identifier, matching a column of
#'     \code{\link{taxa_demo}}.}
#'   \item{subject}{Infant identifier.}
#'   \item{time}{Visit index, 0 to 6.}
#'   \item{month}{Age in months at the visit.}
#'   \item{gender}{Male or Female.}
#'   \item{panc_insuf}{Pancreatic insufficiency, 0 or 1.}
#'   \item{genotype}{Genotype category.}
#'   \item{antibiotics}{Current or prior antibiotics, Y or N.}
#'   \item{breastmilk, formula, table_food}{Feeding at the visit, Y or N.}
#'   \item{weight, length}{Weight in kg and length in cm.}
#' }
#'
#' @source
#' The BONUS cystic-fibrosis infant cohort. Prepared by
#' \code{data-raw/make_taxa_demo.R} in the package sources.
#'
#' @usage data(taxa_demo_meta)
#'
#' @seealso \code{\link{taxa_demo}}, \code{\link{mc_demo_data}}
#'
#' @examples
#' data(taxa_demo_meta)
#'
#' head(taxa_demo_meta[, c("sample", "subject", "time", "month")])
#' table(taxa_demo_meta$gender)
"taxa_demo_meta"
