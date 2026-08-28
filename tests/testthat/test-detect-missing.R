test_that("a complete table reports nothing missing", {
    # every subject-timepoint present and populated
    dat <- data.frame(OTU_ID = c("T1", "T2"))
    for (s in c("S1", "S2")) {
        for (tt in 0:3) dat[[paste0(s, ".", tt)]] <- c(1, 2)
    }

    info <- mc_detect_missing(dat, taxon_col = "OTU_ID")

    expect_equal(nrow(info$missing), 0)
    expect_equal(info$reps, c("S1", "S2"))
    expect_equal(info$times, 0:3)
})

test_that("a column absent from the table is detected", {
    dat <- data.frame(OTU_ID = c("T1", "T2"))
    for (s in c("S1", "S2")) {
        for (tt in 0:3) dat[[paste0(s, ".", tt)]] <- c(1, 2)
    }
    dat[["S1.2"]] <- NULL

    info <- mc_detect_missing(dat, taxon_col = "OTU_ID")

    expect_equal(nrow(info$missing), 1)
    expect_equal(info$missing$rep, "S1")
    expect_equal(info$missing$time, 2)
    expect_equal(info$missing$reason, "absent_column")
    expect_equal(info$missing$col, "S1.2")
})

test_that("a column present but entirely NA is detected", {
    dat <- data.frame(OTU_ID = c("T1", "T2"))
    for (s in c("S1", "S2")) {
        for (tt in 0:3) dat[[paste0(s, ".", tt)]] <- c(1, 2)
    }
    dat[["S2.1"]] <- NA_real_

    info <- mc_detect_missing(dat, taxon_col = "OTU_ID")

    expect_equal(nrow(info$missing), 1)
    expect_equal(info$missing$reason, "all_na")
    expect_equal(info$missing$col, "S2.1")
})

test_that("a partly-NA column is rejected, not silently ignored", {
    # The method imputes whole missing samples. A scattered NA inside an
    # otherwise observed sample is a different problem, and the data has to
    # be corrected before anything is computed from it.
    dat <- data.frame(OTU_ID = c("T1", "T2", "T3"))
    for (s in c("S1", "S2")) {
        for (tt in 0:3) dat[[paste0(s, ".", tt)]] <- c(1, 2, 3)
    }
    dat[["S1.1"]][2] <- NA

    expect_error(
        mc_detect_missing(dat, taxon_col = "OTU_ID"),
        "only partly measured"
    )
})

test_that("an empty taxon does not make every column look partly measured", {
    # A taxon that is NA everywhere accounts for its own missing values. It
    # must not be mistaken for every sample being partly measured.
    dat <- data.frame(OTU_ID = c("T1", "T2", "T3"))
    for (s in c("S1", "S2")) {
        for (tt in 0:3) dat[[paste0(s, ".", tt)]] <- c(1, 2, 3)
    }
    dat[2, grep("\\.", names(dat))] <- NA

    info <- mc_detect_missing(dat, taxon_col = "OTU_ID")
    expect_equal(nrow(info$partial_na), 0)
    expect_equal(nrow(info$missing), 0)
})

test_that("an all-NA column of logical type is accepted", {
    # read.csv() returns an all-NA column as logical, not numeric. That is a
    # normal way for a missing sample to arrive and must not be rejected.
    dat <- data.frame(OTU_ID = c("T1", "T2"))
    for (s in c("S1", "S2")) {
        for (tt in 0:3) dat[[paste0(s, ".", tt)]] <- c(1, 2)
    }
    dat[["S2.1"]] <- NA # logical, not NA_real_
    expect_type(dat[["S2.1"]], "logical")

    info <- mc_detect_missing(dat, taxon_col = "OTU_ID")

    expect_equal(nrow(info$missing), 1)
    expect_equal(info$missing$reason, "all_na")
})

test_that("the bundled CSV round-trips through read.csv", {
    path <- system.file("extdata", "taxa_demo.csv", package = "MicrobiomeCurves")
    skip_if(path == "", "extdata not installed")

    dat <- read.csv(path, check.names = FALSE)
    info <- mc_detect_missing(dat, taxon_col = "OTU_ID")

    expect_equal(nrow(info$missing), 3)
    expect_setequal(info$missing$col, c("S02.1", "S04.2", "S07.4"))
})

test_that("a non-numeric column is rejected", {
    dat <- data.frame(OTU_ID = c("T1", "T2"))
    for (s in c("S1", "S2")) {
        for (tt in 0:3) dat[[paste0(s, ".", tt)]] <- c(1, 2)
    }
    dat[["S1.1"]] <- c("a", "b")

    expect_error(
        mc_detect_missing(dat, taxon_col = "OTU_ID"),
        "not numeric"
    )
})

test_that("columns that do not parse are reported, not silently used", {
    # Real sorted tables carry metadata columns such as `cluster`.
    dat <- data.frame(OTU_ID = c("T1", "T2"))
    for (s in c("S1", "S2")) {
        for (tt in 0:3) dat[[paste0(s, ".", tt)]] <- c(1, 2)
    }
    dat[["cluster"]] <- c(1, 2)

    info <- mc_detect_missing(dat, taxon_col = "OTU_ID")

    expect_equal(info$unparsed, "cluster")
    expect_false("cluster" %in% info$col_map$col)
})

test_that("observed counts per subject are correct", {
    dat <- data.frame(OTU_ID = c("T1", "T2"))
    for (s in c("S1", "S2")) {
        for (tt in 0:3) dat[[paste0(s, ".", tt)]] <- c(1, 2)
    }
    dat[["S1.2"]] <- NULL # absent
    dat[["S1.3"]] <- NA_real_ # all NA

    info <- mc_detect_missing(dat, taxon_col = "OTU_ID")

    expect_equal(info$observed$n_observed[info$observed$rep == "S1"], 2)
    expect_equal(info$observed$n_observed[info$observed$rep == "S2"], 4)
})

test_that("taxon_col and layout are validated", {
    expect_error(mc_detect_missing(list(a = 1)), "data.frame")
    expect_error(
        mc_detect_missing(data.frame(x = 1), taxon_col = "nope"),
        "not found"
    )
    expect_error(
        mc_detect_missing(data.frame(OTU_ID = "T1", junk = 1)),
        "No '<subject>[.]<time>' columns detected"
    )
})
