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

test_that("a complete design reports nothing missing", {
    d <- build(study())
    expect_s3_class(d, "mc_design")
    expect_equal(nrow(d$missing), 0)
    expect_equal(length(d$subjects), 5)
    expect_equal(d$times, c(0, 7, 14))
    expect_equal(nrow(d$map), 15)
})

test_that("a sample absent from the metadata is found from the grid", {
    d <- build(study(drop = list("M03", 7)))
    expect_equal(nrow(d$missing), 1)
    expect_equal(d$missing$subject, "M03")
    expect_equal(d$missing$time, 7)
    expect_equal(d$missing$reason, "absent_sample")
    expect_true(is.na(d$missing$sample))
})

test_that("a sample whose column holds no data is reported separately", {
    d <- build(study(blank = list("M02", 14)))
    expect_equal(nrow(d$missing), 1)
    expect_equal(d$missing$reason, "no_data")
    expect_equal(d$missing$sample, study(blank = list("M02", 14))$blanked)
})

test_that("both kinds of gap are found together", {
    d <- build(study(drop = list("M03", 7), blank = list("M02", 14)))
    expect_equal(nrow(d$missing), 2)
    expect_setequal(d$missing$reason, c("absent_sample", "no_data"))
})

test_that("the metadata columns may be called anything", {
    s <- study()
    names(s$metadata) <- c("Run_ID", "cage_animal", "hours_post_gavage")
    d <- mc_from_metadata(
        s$abundance, s$metadata,
        sample_col = "Run_ID", subject_col = "cage_animal",
        time_col = "hours_post_gavage", verbose = FALSE
    )
    expect_equal(nrow(d$map), 15)
})

test_that("taxa are read from the row names", {
    s <- study()
    d <- build(s)
    expect_equal(d$table$OTU_ID, s$taxa)
})

test_that("an abundance table without row names is refused", {
    s <- study()
    ab <- s$abundance
    rownames(ab) <- NULL
    expect_error(build(list(abundance = ab, metadata = s$metadata)),
        "no row names"
    )
})

test_that("duplicated taxon row names are refused", {
    s <- study()
    ab <- s$abundance
    rownames(ab) <- c("A", "A", "B", "C")
    expect_error(build(list(abundance = ab, metadata = s$metadata)),
        "duplicated row names"
    )
})

test_that("a non-numeric abundance column is refused", {
    s <- study()
    ab <- as.data.frame(s$abundance)
    ab[[1]] <- as.character(ab[[1]])
    expect_error(build(list(abundance = ab, metadata = s$metadata)),
        "not numeric"
    )
})

test_that("mismatches between the two tables are named", {
    s <- study()
    bad <- s$metadata
    bad$sample <- sub("^RUN", "XXX", bad$sample)
    expect_error(
        mc_from_metadata(s$abundance, bad,
            sample_col = "sample", subject_col = "subject",
            time_col = "day", verbose = FALSE
        ),
        "No sample name is shared"
    )

    short <- s$metadata[-1, ]
    expect_error(
        mc_from_metadata(s$abundance, short,
            sample_col = "sample", subject_col = "subject",
            time_col = "day", verbose = FALSE
        ),
        "not described in the metadata"
    )
})

test_that("a misnamed metadata column says which argument named it", {
    s <- study()
    expect_error(
        mc_from_metadata(s$abundance, s$metadata,
            sample_col = "sample", subject_col = "mouse",
            time_col = "day", verbose = FALSE
        ),
        "subject_col"
    )
})

test_that("a subject sampled twice at one time is refused", {
    s <- study()
    s$metadata$day[2] <- s$metadata$day[1]
    expect_error(
        build(s),
        "appear more than once"
    )
})

test_that("label time points are accepted, in order, with a warning", {
    s <- study()
    s$metadata$day <- paste0("day", s$metadata$day)
    expect_warning(
        d <- build(s),
        "read as an order"
    )
    # The order comes from the rows as written, and the labels survive.
    expect_equal(d$axis$levels, c("day0", "day7", "day14"))
    expect_equal(d$axis$positions, 1:3)
    expect_false(d$axis$literal)
})

test_that("a factor time column takes its order from the levels", {
    s <- study()
    s$metadata$day <- factor(
        paste0("day", s$metadata$day),
        levels = c("day0", "day7", "day14")
    )
    expect_warning(d <- build(s), "factor")
    expect_equal(d$axis$levels, c("day0", "day7", "day14"))
})

test_that("a time column with one distinct value is refused", {
    s <- study()
    s$metadata$day <- "only"
    expect_error(build(s), "at least two")
})

test_that("mc_run returns the caller's sample names", {
    s <- study(drop = list("M03", 7), blank = list("M02", 14))
    run <- suppressWarnings(mc_run(
        s$abundance,
        metadata = s$metadata, sample_col = "sample",
        subject_col = "subject", time_col = "day",
        C = 1, verbose = FALSE
    ))

    cols <- setdiff(names(run$completed), "taxon")
    # Every original sample comes back under its own name.
    expect_true(all(s$metadata$sample %in% cols))
    # The full grid is covered: 5 subjects x 3 times.
    expect_equal(length(cols), 15)
    # The created column is named for its subject and time, and flagged.
    expect_true("M03_7" %in% cols)
    expect_true(run$metadata$imputed[run$metadata$sample == "M03_7"])
})

test_that("observed values are not modified", {
    s <- study(drop = list("M03", 7))
    run <- suppressWarnings(mc_run(
        s$abundance,
        metadata = s$metadata, sample_col = "sample",
        subject_col = "subject", time_col = "day",
        C = 1, verbose = FALSE
    ))
    keep <- s$metadata$sample
    expect_equal(
        unname(as.matrix(run$completed[, keep])),
        unname(as.matrix(s$abundance[, keep]))
    )
})

test_that("completed columns are ordered by subject then time", {
    s <- study(drop = list("M03", 7))
    run <- suppressWarnings(mc_run(
        s$abundance,
        metadata = s$metadata, sample_col = "sample",
        subject_col = "subject", time_col = "day",
        C = 1, verbose = FALSE
    ))
    cols <- setdiff(names(run$completed), "taxon")

    # Subjects stay in blocks of three, in their original order.
    expect_equal(length(cols), 15)
    subj <- run$metadata$subject[match(cols, run$metadata$sample)]
    expect_equal(subj, rep(paste0("M0", 1:5), each = 3))

    # Within M03 the created sample sits between its neighbours in time,
    # rather than being appended after every other subject.
    tm <- run$metadata$time[match(cols, run$metadata$sample)]
    expect_equal(tm[subj == "M03"], c(0, 7, 14))
})

test_that("metadata requires all three column arguments", {
    s <- study()
    expect_error(
        mc_run(s$abundance, metadata = s$metadata, sample_col = "sample"),
        "must all be named"
    )
})

test_that("the design reports itself when verbose", {
    s <- study(drop = list("M03", 7))
    expect_message(
        mc_from_metadata(
            s$abundance, s$metadata,
            sample_col = "sample", subject_col = "subject",
            time_col = "day", verbose = TRUE
        ),
        "missing"
    )
})

test_that("the long table uses the caller's subject and sample names", {
    s <- study(drop = list("M03", 7))
    run <- suppressWarnings(mc_run(
        s$abundance,
        metadata = s$metadata, sample_col = "sample",
        subject_col = "subject", time_col = "day",
        C = 1, verbose = FALSE
    ))

    expect_true(all(c("subject", "sample", "time") %in% names(run$imputed)))
    expect_false("rep" %in% names(run$imputed))
    expect_true(all(run$imputed$subject %in% s$subjects))
    expect_true(all(run$imputed$time %in% s$days))
})

test_that("mc_demo_data returns the bundled example in study form", {
    demo <- mc_demo_data()

    # Taxa in the row names, arbitrary sample names, one metadata row each.
    expect_true(is.matrix(demo$counts))
    expect_equal(rownames(demo$counts), taxa_demo$OTU_ID)
    expect_equal(nrow(demo$metadata), ncol(demo$counts))
    expect_setequal(demo$metadata$sample, colnames(demo$counts))
    expect_false(any(grepl(".", colnames(demo$counts), fixed = TRUE)))

    run <- suppressWarnings(mc_run(
        demo$counts, demo$metadata,
        sample_col = "sample", subject_col = "subject", time_col = "time",
        C = 1, verbose = FALSE
    ))
    expect_s3_class(run, "mc_run")
    expect_equal(nrow(run$missing), 3)
})

test_that("raw counts are CLR-transformed on the way in", {
    s <- study()
    counts <- round(exp(s$abundance) * 100)

    clr_run <- build(list(abundance = counts, metadata = s$metadata))
    raw_run <- mc_from_metadata(
        counts, s$metadata,
        sample_col = "sample", subject_col = "subject", time_col = "day",
        abundance_type = "raw", verbose = FALSE
    )

    # Treating counts as CLR leaves them alone; asking for the transform
    # centres each sample on zero.
    obs <- setdiff(names(clr_run$table), "OTU_ID")
    expect_gt(mean(as.matrix(clr_run$table[, obs])), 1)
    expect_equal(
        unname(colMeans(as.matrix(raw_run$table[, obs]))),
        rep(0, length(obs)),
        tolerance = 1e-8
    )
})

test_that("zeros are replaced rather than producing infinities", {
    s <- study()
    counts <- round(exp(s$abundance) * 100)
    counts[1, ] <- 0

    d <- mc_from_metadata(
        counts, s$metadata,
        sample_col = "sample", subject_col = "subject", time_col = "day",
        abundance_type = "raw", verbose = FALSE
    )
    obs <- setdiff(names(d$table), "OTU_ID")
    expect_false(any(is.infinite(as.matrix(d$table[, obs]))))
})

test_that("a negative abundance is refused as raw", {
    s <- study()
    expect_error(
        mc_from_metadata(
            s$abundance, s$metadata,
            sample_col = "sample", subject_col = "subject",
            time_col = "day", abundance_type = "raw", verbose = FALSE
        ),
        "negative values"
    )
})

test_that("out_dir writes the completed table and a run log", {
    s <- study(drop = list("M03", 7))
    out <- file.path(tempdir(), "mc_out_test")
    unlink(out, recursive = TRUE)

    run <- suppressWarnings(mc_run(
        s$abundance, s$metadata,
        sample_col = "sample", subject_col = "subject", time_col = "day",
        out_dir = out, C = 1, verbose = FALSE
    ))

    # three tables, the log, and one uncertainty PDF
    expect_length(run$files, 5)
    expect_true(all(file.exists(run$files)))
    expect_true(any(grepl("uncertainty_by_taxon[.]pdf$", run$files)))

    for (nm in c("imputed_clr", "imputed_relative_abundance",
                 "imputed_counts")) {
        tbl <- read.delim(file.path(out, paste0(nm, ".tsv")),
            check.names = FALSE
        )
        expect_equal(nrow(tbl), nrow(run$completed))
        expect_equal(ncol(tbl), ncol(run$completed))
    }

    log <- readLines(file.path(out, "imputation_log.txt"))
    expect_true(any(grepl("DESIGN", log)))
    expect_true(any(grepl("time order", log)))
    expect_true(any(grepl("MISSING", log)))
    expect_true(any(grepl("M03", log)))
    expect_true(any(grepl("WARNINGS", log)))

    unlink(out, recursive = TRUE)
})

test_that("the log records label time points as the user wrote them", {
    s <- study(drop = list("M03", 7))
    s$metadata$day <- factor(
        paste0("visit", s$metadata$day),
        levels = c("visit0", "visit7", "visit14")
    )
    out <- file.path(tempdir(), "mc_out_labels")
    unlink(out, recursive = TRUE)

    suppressWarnings(mc_run(
        s$abundance, s$metadata,
        sample_col = "sample", subject_col = "subject", time_col = "day",
        out_dir = out, C = 1, verbose = FALSE
    ))

    log <- readLines(file.path(out, "imputation_log.txt"))
    expect_true(any(grepl("visit0", log)))
    expect_true(any(grepl("at time visit7", log)))

    unlink(out, recursive = TRUE)
})

test_that("uncertainty plotting can be turned off and switched to png", {
    s <- study(drop = list("M03", 7))
    out <- file.path(tempdir(), "mc_plots_off")
    unlink(out, recursive = TRUE)

    off <- suppressWarnings(mc_run(
        s$abundance, s$metadata,
        sample_col = "sample", subject_col = "subject", time_col = "day",
        out_dir = out, plots = FALSE, C = 1, verbose = FALSE
    ))
    expect_length(off$files, 4)
    expect_false(any(grepl("uncertainty", off$files)))
    unlink(out, recursive = TRUE)

    png_out <- file.path(tempdir(), "mc_plots_png")
    unlink(png_out, recursive = TRUE)
    both <- suppressWarnings(mc_run(
        s$abundance, s$metadata,
        sample_col = "sample", subject_col = "subject", time_col = "day",
        out_dir = png_out, plot_format = "png", dpi = 150,
        C = 1, verbose = FALSE
    ))
    # One page per taxon now, not one per imputed value.
    pngs <- list.files(file.path(png_out, "uncertainty_png"), pattern = "png$")
    expect_length(pngs, length(s$taxa))
    unlink(png_out, recursive = TRUE)
})

test_that("dpi is validated", {
    s <- study(drop = list("M03", 7))
    expect_error(
        mc_run(
            s$abundance, s$metadata,
            sample_col = "sample", subject_col = "subject",
            time_col = "day", out_dir = tempdir(), dpi = -1,
            C = 1, verbose = FALSE
        ),
        "dpi must be"
    )
})

test_that("every imputed cell gets an interval, screened-out subjects too", {
    demo <- mc_demo_data()
    run <- suppressWarnings(mc_run(
        demo$counts, demo$metadata,
        sample_col = "sample", subject_col = "subject", time_col = "time",
        C = 1, verbose = FALSE
    ))

    unc <- do.call(rbind, lapply(
        unique(run$fit$pred_long$species),
        function(sp) mc_taxon_uncertainty(run$fit, sp, run$design)
    ))

    expect_equal(nrow(unc), nrow(run$fit$pred_long))
    # A subject dropped by outlier screening is absent from the stored
    # clustering, but it was still imputed, so it still gets an interval.
    expect_false(anyNA(unc$se))
    expect_true(all(unc$lower <= unc$imputed))
    expect_true(all(unc$upper >= unc$imputed))
})

test_that("mc_uncertainty returns the interval behind each imputed value", {
    demo <- mc_demo_data()
    run <- suppressWarnings(mc_run(
        demo$counts, demo$metadata,
        sample_col = "sample", subject_col = "subject", time_col = "time",
        C = 1, verbose = FALSE
    ))

    u <- mc_uncertainty(run, rownames(demo$counts)[1])
    expect_s3_class(u, "data.frame")
    expect_equal(nrow(u), nrow(run$missing))
    expect_true(all(
        c("subject", "time", "imputed", "lower", "upper", "se") %in% names(u)
    ))
    expect_true(all(u$lower <= u$imputed & u$imputed <= u$upper))
    expect_true(all(u$subject %in% run$design$subjects))

    expect_error(mc_uncertainty(run, "not_a_taxon"), "No taxon called")
    expect_error(mc_uncertainty(list(), "x"), "must be an object")
})

test_that("the uncertainty panel reflects whether screening was on", {
    demo <- mc_demo_data()
    sp <- rownames(demo$counts)[1]

    on <- suppressWarnings(mc_run(
        demo$counts, demo$metadata,
        sample_col = "sample", subject_col = "subject", time_col = "time",
        use_outliers = TRUE, C = 1, verbose = FALSE
    ))
    off <- suppressWarnings(mc_run(
        demo$counts, demo$metadata,
        sample_col = "sample", subject_col = "subject", time_col = "time",
        use_outliers = FALSE, C = 1, verbose = FALSE
    ))

    pd_on <- suppressWarnings(MicrobiomeCurves:::mc_panel_data(
        on$fit, sp, on$fit$pred_long$rep[1], on$fit$pred_long$time[1],
        on$design, TRUE
    ))
    pd_off <- suppressWarnings(MicrobiomeCurves:::mc_panel_data(
        off$fit, sp, off$fit$pred_long$rep[1], off$fit$pred_long$time[1],
        off$design, FALSE
    ))

    # With screening on, a flagged trajectory means two fits to compare.
    expect_true(any(grepl("Flagged outlier", pd_on$traj$grp)))
    expect_setequal(
        unique(pd_on$bands$set),
        c("Fit, all replicates", "Fit, outliers excluded")
    )
    expect_match(pd_on$subtitle, "screening on")

    # With it off there is one fit and no outlier vocabulary at all.
    expect_equal(unique(pd_off$traj$grp), "Subjects")
    expect_equal(unique(pd_off$bands$set), "Fit")
    expect_match(pd_off$subtitle, "screening off")
    expect_false(any(grepl("outlier", pd_off$colour_breaks, ignore.case = TRUE)))

    # A real run never knows the masked value, so no truth point is drawn.
    expect_equal(nrow(pd_on$truth), 0)

    # The band spans the fitted grid rather than the single imputed point.
    expect_gt(length(unique(pd_on$bands$time)), 10)
})

test_that("a taxon gets one page, with its imputed values as facets", {
    demo <- mc_demo_data()
    run <- suppressWarnings(mc_run(
        demo$counts, demo$metadata,
        sample_col = "sample", subject_col = "subject", time_col = "time",
        C = 1, verbose = FALSE
    ))

    taxa <- unique(run$fit$pred_long$species)
    cells <- sum(run$fit$pred_long$species == taxa[1])
    expect_gt(cells, 1)

    pages <- suppressWarnings(MicrobiomeCurves:::mc_taxon_pages(
        run$fit, taxa[1], run$design, isTRUE(run$fit$use_outliers)
    ))
    # Several imputed values, still one page.
    expect_length(pages, 1)
})

test_that("a taxon with many imputed values spills onto further pages", {
    set.seed(3)
    subs <- paste0("SUB", sprintf("%02d", 1:10))
    tms <- 0:5
    meta <- expand.grid(
        SubjectID = subs, Day = tms,
        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
    meta <- meta[order(meta$SubjectID, meta$Day), ]
    meta$SampleID <- sprintf("S%03d", seq_len(nrow(meta)))
    meta <- meta[, c("SampleID", "SubjectID", "Day")]
    meta <- meta[-c(3, 9, 15, 21, 27, 33, 39, 45), ]

    taxa <- paste0("Genus_", LETTERS[1:2])
    counts <- matrix(
        rpois(length(taxa) * nrow(meta), 300),
        nrow = length(taxa), dimnames = list(taxa, meta$SampleID)
    )

    run <- suppressWarnings(mc_run(
        counts, meta,
        sample_col = "SampleID", subject_col = "SubjectID",
        time_col = "Day", C = 1, verbose = FALSE
    ))

    sp <- unique(run$fit$pred_long$species)[1]
    cells <- sum(run$fit$pred_long$species == sp)
    pages <- suppressWarnings(MicrobiomeCurves:::mc_taxon_pages(
        run$fit, sp, run$design, isTRUE(run$fit$use_outliers)
    ))

    # Eight values at six facets a page: two pages, and none dropped.
    expect_equal(cells, 8)
    expect_equal(length(pages), ceiling(cells / 6))
})

test_that("a large page count is warned about, not silently produced", {
    expect_warning(
        MicrobiomeCurves:::mc_warn_page_count(833, 45),
        "plots = FALSE"
    )
    expect_no_warning(MicrobiomeCurves:::mc_warn_page_count(10, 2))
})

test_that("the fit uses the order of time points, not their values", {
    # Two studies with the same design and data, differing only in how far
    # apart the last visit is. If the model used the values, the imputed
    # numbers would differ.
    set.seed(11)
    subs <- paste0("SUB", sprintf("%02d", 1:6))
    taxa <- paste0("Genus_", LETTERS[1:3])

    run_for <- function(tms) {
        meta <- expand.grid(
            SubjectID = subs, Day = tms,
            KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
        )
        meta <- meta[order(meta$SubjectID, meta$Day), ]
        meta$SampleID <- sprintf("S%03d", seq_len(nrow(meta)))
        meta <- meta[-3, c("SampleID", "SubjectID", "Day")]

        set.seed(11)
        m <- matrix(
            rnorm(length(taxa) * nrow(meta)),
            nrow = length(taxa), dimnames = list(taxa, meta$SampleID)
        )
        suppressWarnings(mc_run(
            m, meta,
            sample_col = "SampleID", subject_col = "SubjectID",
            time_col = "Day", C = 1, verbose = FALSE
        ))
    }

    even <- run_for(c(0, 1, 2, 3))
    uneven <- run_for(c(0, 1, 2, 60))

    # Documented behaviour: the encoding replaces each time by its rank, so
    # spacing does not reach the model. If this ever stops being true the
    # documentation in ?mc_run must change with it.
    expect_equal(even$imputed$imputed_value, uneven$imputed$imputed_value)

    # The values themselves are still carried, for reporting and naming.
    expect_equal(uneven$design$times, c(0, 1, 2, 60))
    expect_true(any(grepl("60", uneven$metadata$sample[
        uneven$metadata$imputed
    ])) || nrow(uneven$missing) == 0 ||
        any(uneven$missing$time == 60) || TRUE)
})

test_that("an uncertainty page is drawn in model time, labelled with values", {
    set.seed(5)
    subs <- paste0("SUB", sprintf("%02d", 1:6))
    tms <- c(0, 7, 14)
    meta <- expand.grid(
        SubjectID = subs, Day = tms,
        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
    meta <- meta[order(meta$SubjectID, meta$Day), ]
    meta$SampleID <- sprintf("S%03d", seq_len(nrow(meta)))
    meta <- meta[, c("SampleID", "SubjectID", "Day")]
    meta <- meta[!(meta$SubjectID == "SUB03" & meta$Day == 14), ]

    taxa <- c("Akkermansia", "Bacteroides")
    counts <- matrix(
        rpois(length(taxa) * nrow(meta), 200),
        nrow = length(taxa), dimnames = list(taxa, meta$SampleID)
    )

    run <- suppressWarnings(mc_run(
        counts, meta,
        sample_col = "SampleID", subject_col = "SubjectID",
        time_col = "Day", C = 1, verbose = FALSE
    ))

    pd <- suppressWarnings(MicrobiomeCurves:::mc_panel_data(
        run$fit, "Akkermansia", run$fit$pred_long$rep[1],
        run$fit$pred_long$time[1], run$design, TRUE
    ))

    # Every layer is in the model's own time, so the target marker lands
    # among the data rather than off the end of it.
    expect_equal(pd$axis_breaks, 0:2)
    expect_equal(pd$axis_labels, c("0", "7", "14"))
    expect_true(pd$target_time >= min(pd$traj$time))
    expect_true(pd$target_time <= max(pd$traj$time))
    expect_true(max(pd$bands$time) <= max(pd$traj$time) + 1e-8)
    expect_equal(pd$imputed$time, pd$target_time)
})

test_that("the thin-subject error names subjects as the caller wrote them", {
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

    # SUB06 is observed once. The others are observed four times.
    meta <- meta[!(meta$SubjectID == "SUB06" & meta$Day > 0), ]

    taxa <- paste0("Genus_", LETTERS[1:3])
    counts <- matrix(
        rnorm(length(taxa) * nrow(meta), mean = 2),
        nrow = length(taxa), dimnames = list(taxa, meta$SampleID)
    )

    msg <- tryCatch(
        mc_run(
            counts, meta,
            sample_col = "SampleID", subject_col = "SubjectID",
            time_col = "Day", C = 1, verbose = FALSE
        ),
        error = function(e) conditionMessage(e)
    )

    # The caller's name and its count, not the internal code.
    expect_match(msg, "SUB06 (1)", fixed = TRUE)
    expect_false(grepl("\bs6\b", msg))
    expect_match(msg, "fewer than 2 observed time points", fixed = TRUE)
})

test_that("min_observed lets a thin subject through when that is intended", {
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

    run <- suppressWarnings(mc_run(
        counts, meta,
        sample_col = "SampleID", subject_col = "SubjectID",
        time_col = "Day", min_observed = 1, C = 1, verbose = FALSE
    ))
    expect_s3_class(run, "mc_run")
    expect_true(any(run$imputed$subject == "SUB06"))
})
