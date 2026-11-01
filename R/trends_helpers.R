############################################
# Internals for mc_trends() and mc_compare().
#
# Both functions fit the same model and differ in one term, so everything
# except that term lives here and is shared. The model is a generalised
# additive model: abundance against a smooth of time, with a penalised
# random effect for subject so that repeated measures are not treated as
# independent samples.
#
# Nothing here touches the imputation path. A completed table from mc_run()
# is one valid input, an abundance table with gaps is another: a GAM drops
# incomplete rows and fits on what is left, so the caller chooses whether
# imputed cells take part by choosing which table to pass.
############################################


#' Column names this file reserves inside the per-taxon data frame
#'
#' @description
#' The formula is built as text, so the columns it names have to be called
#' the same thing for every taxon whatever the caller's own columns are
#' called. These four names are therefore taken, and a covariate using one
#' of them stops the call rather than silently shadowing it.
#'
#' @keywords internal
#' @noRd
mc_trend_reserved <- c("value", "time", "subject", "mc_grp")


#' Build the per-sample design for a trend model
#'
#' @param metadata The caller's metadata.
#' @param sample_col,subject_col,time_col Column names within `metadata`.
#' @param group_col Grouping column, or `NULL`.
#' @param covariates Character vector of columns to adjust for, or `NULL`.
#'
#' @return A data.frame with `sample`, `subject`, `time`, the covariates
#'   under their own names, and `mc_grp` when `group_col` was given.
#'
#' @keywords internal
#' @noRd
mc_trend_design <- function(metadata, sample_col, subject_col, time_col,
                            group_col = NULL, covariates = NULL) {
    if (!is.data.frame(metadata)) {
        stop("metadata must be a data.frame.", call. = FALSE)
    }

    need <- c(sample_col, subject_col, time_col, group_col, covariates)
    gone <- setdiff(need, names(metadata))
    if (length(gone) > 0) {
        stop(
            "metadata has no column(s) ", mc_fmt_some(gone), ".",
            call. = FALSE
        )
    }

    bad <- intersect(covariates, mc_trend_reserved)
    if (length(bad) > 0) {
        stop(
            "The covariate name(s) ", mc_fmt_some(bad), " are used ",
            "internally by the model formula. Rename the column(s) before ",
            "calling.",
            call. = FALSE
        )
    }

    time <- metadata[[time_col]]
    # A smooth needs a continuous axis. Coercing a label such as "week_4"
    # would invent an ordering and equal spacing that the study may not have,
    # so it is refused rather than guessed.
    if (!is.numeric(time)) {
        stop(
            "time_col '", time_col, "' is ", class(time)[1], ", not numeric. ",
            "s(time) needs a continuous axis: convert the column to the ",
            "elapsed time it stands for (days, weeks) before calling.",
            call. = FALSE
        )
    }

    d <- data.frame(
        sample = as.character(metadata[[sample_col]]),
        subject = as.character(metadata[[subject_col]]),
        time = as.numeric(time),
        stringsAsFactors = FALSE
    )

    if (!is.null(group_col)) {
        d$mc_grp <- as.character(metadata[[group_col]])
    }
    for (cv in covariates) {
        d[[cv]] <- metadata[[cv]]
    }

    keep <- !is.na(d$time) & !is.na(d$subject) & nzchar(d$subject)
    if (!is.null(group_col)) keep <- keep & !is.na(d$mc_grp)
    d <- d[keep, , drop = FALSE]

    if (anyDuplicated(d$sample) > 0) {
        dup <- unique(d$sample[duplicated(d$sample)])
        stop(
            "metadata has more than one row for sample(s) ",
            mc_fmt_some(dup), ".",
            call. = FALSE
        )
    }
    if (nrow(d) == 0) {
        stop("No usable rows in metadata after dropping NA.", call. = FALSE)
    }
    d
}


#' Pull the abundance table apart into a matrix and taxon labels
#'
#' @param abundance A data.frame with a taxon column, or a matrix with
#'   row names.
#' @param taxon_col Name of the taxon column, or `NULL` to detect it.
#' @param taxa Character vector of taxa to keep, or `NULL` for all.
#'
#' @return A list of `mat` (taxa by samples, numeric) and `taxa`.
#'
#' @keywords internal
#' @noRd
mc_trend_matrix <- function(abundance, taxon_col = NULL, taxa = NULL) {
    if (is.matrix(abundance)) {
        if (is.null(rownames(abundance))) {
            stop(
                "A matrix abundance needs row names to identify the taxa.",
                call. = FALSE
            )
        }
        mat <- abundance
        labs <- rownames(abundance)
    } else if (is.data.frame(abundance)) {
        if (is.null(taxon_col)) {
            # The completed table from mc_run() calls it "taxon"; a raw table
            # usually calls it "OTU_ID" or similar. Either way it is the one
            # column that does not hold numbers.
            chr <- which(!vapply(abundance, is.numeric, logical(1)))
            if (length(chr) == 0) {
                stop(
                    "abundance has no non-numeric column to read taxon ",
                    "names from. Name it with taxon_col, or supply a matrix ",
                    "with row names.",
                    call. = FALSE
                )
            }
            taxon_col <- names(abundance)[chr[1]]
        }
        if (!taxon_col %in% names(abundance)) {
            stop(
                "abundance has no column '", taxon_col, "'.",
                call. = FALSE
            )
        }
        labs <- as.character(abundance[[taxon_col]])
        mat <- as.matrix(abundance[, setdiff(names(abundance), taxon_col),
            drop = FALSE
        ])
    } else {
        stop(
            "abundance must be a data.frame or a matrix.",
            call. = FALSE
        )
    }

    storage.mode(mat) <- "numeric"
    rownames(mat) <- labs

    if (!is.null(taxa)) {
        gone <- setdiff(taxa, labs)
        if (length(gone) > 0) {
            stop(
                "abundance has no row(s) for taxa ", mc_fmt_some(gone), ".",
                call. = FALSE
            )
        }
        mat <- mat[match(taxa, labs), , drop = FALSE]
        labs <- taxa
    }
    if (anyDuplicated(labs) > 0) {
        stop(
            "abundance has repeated taxon name(s): ",
            mc_fmt_some(unique(labs[duplicated(labs)])), ".",
            call. = FALSE
        )
    }
    list(mat = mat, taxa = labs)
}


#' Keep only the samples that the table and the metadata agree on
#'
#' @param design Per-sample design from [mc_trend_design()].
#' @param mat Abundance matrix from [mc_trend_matrix()].
#' @param verbose Logical. Whether to report what was dropped.
#'
#' @return A list of the aligned `design` and `mat`.
#'
#' @keywords internal
#' @noRd
mc_trend_align <- function(design, mat, verbose = TRUE) {
    shared <- intersect(design$sample, colnames(mat))
    if (length(shared) == 0) {
        stop(
            "No sample name is in both abundance and metadata. The ",
            "abundance column names and metadata[[sample_col]] must agree.",
            call. = FALSE
        )
    }

    only_meta <- setdiff(design$sample, colnames(mat))
    only_tab <- setdiff(colnames(mat), design$sample)
    if (isTRUE(verbose) && (length(only_meta) > 0 || length(only_tab) > 0)) {
        message(
            "  using ", length(shared), " shared sample(s); dropped ",
            length(only_meta), " in metadata only, ",
            length(only_tab), " in abundance only"
        )
    }

    design <- design[match(shared, design$sample), , drop = FALSE]
    rownames(design) <- NULL
    list(design = design, mat = mat[, shared, drop = FALSE])
}


#' Resolve the basis size
#'
#' @description
#' mgcv defaults to `k = 10`, which needs ten distinct values of the
#' covariate and fails outright below that. Most longitudinal microbiome
#' designs have fewer time points than that, so the default here is read
#' from the design instead.
#'
#' @param times Numeric vector of observed times.
#' @param k The caller's `k`, or `NULL`.
#'
#' @return An integer basis size.
#'
#' @keywords internal
#' @noRd
mc_trend_k <- function(times, k = NULL) {
    n <- length(unique(times[!is.na(times)]))
    if (n < 3L) {
        stop(
            "Only ", n, " distinct time point(s). A smooth of time needs at ",
            "least 3.",
            call. = FALSE
        )
    }
    if (is.null(k)) {
        return(as.integer(min(n - 1L, 10L)))
    }
    if (length(k) != 1L || !is.finite(k) || k < 3L) {
        stop("k must be a single number of at least 3.", call. = FALSE)
    }
    k <- as.integer(k)
    if (k > n) {
        stop(
            "k = ", k, " exceeds the ", n, " distinct time point(s) in the ",
            "data. Lower k, or leave it NULL to have it chosen.",
            call. = FALSE
        )
    }
    k
}


#' The right-hand side terms of a trend model
#'
#' @param k Basis size.
#' @param subject_effect One of `"intercept"`, `"slope"`, `"smooth"`,
#'   `"none"`.
#' @param covariates Character vector, or `NULL`.
#' @param group One of `"none"`, `"ordered"` or `"global"`.
#' @param group_smooth Logical. Whether the group gets a smooth of its own.
#'   `FALSE` builds the null model for the nested comparison.
#'
#' @return A character vector of terms.
#'
#' @keywords internal
#' @noRd
mc_trend_terms <- function(k, subject_effect, covariates = NULL,
                           group = "none", group_smooth = TRUE) {
    tt <- sprintf("s(time, k = %d)", k)

    # The group enters in one of two parameterisations. With an ordered
    # factor the smooth is a difference from the reference level, and the
    # factor's own main effect has to be present for that to be identifiable.
    # With an unordered factor the by-smooths are centred, so s(time) becomes
    # the trend over all groups and each smooth is a departure from it; the
    # factor is then carried parametrically, because three or four arms are
    # too few levels to estimate a random-effect variance from.
    if (identical(group, "ordered")) {
        tt <- c("mc_grp", tt)
        if (isTRUE(group_smooth)) {
            tt <- c(tt, sprintf("s(time, by = mc_grp, k = %d)", k))
        }
    } else if (identical(group, "global")) {
        tt <- c(tt, "mc_grp")
        if (isTRUE(group_smooth)) {
            tt <- c(tt, sprintf("s(time, by = mc_grp, k = %d)", k))
        }
    }

    if (length(covariates) > 0) {
        tt <- c(tt, covariates)
    }

    tt <- c(tt, switch(subject_effect,
        intercept = "s(subject, bs = \"re\")",
        slope = c("s(subject, bs = \"re\")", "s(subject, time, bs = \"re\")"),
        smooth = sprintf("s(time, subject, bs = \"fs\", k = %d, m = 1)", k),
        none = character(0),
        stop("Unknown subject_effect '", subject_effect, "'.", call. = FALSE)
    ))

    tt
}


#' Assemble a formula from terms
#'
#' @param terms Character vector from [mc_trend_terms()].
#'
#' @return A formula with `value` as the response.
#'
#' @keywords internal
#' @noRd
mc_trend_formula <- function(terms) {
    stats::as.formula(paste("value ~", paste(terms, collapse = " + ")))
}


#' Put the group column into the factor form the parameterisation needs
#'
#' @param d The per-taxon data frame.
#' @param levels_in_order Character vector of group levels, reference first.
#' @param group One of `"ordered"` or `"global"`.
#'
#' @return `d`, with `mc_grp` as a factor.
#'
#' @keywords internal
#' @noRd
mc_trend_group_factor <- function(d, levels_in_order, group) {
    d$mc_grp <- factor(d$mc_grp, levels = levels_in_order)
    if (identical(group, "ordered")) {
        d$mc_grp <- factor(d$mc_grp, levels = levels_in_order, ordered = TRUE)
    }
    d
}


#' Resolve which level is the reference, or that there is none
#'
#' @param levs Character vector of the group's levels, sorted.
#' @param reference The caller's `reference`, or `NULL`.
#'
#' @return A list of `group` (`"ordered"` or `"global"`), `reference` (a
#'   level name or `"global"`), `levels` (reference first), and `set_by`.
#'
#' @keywords internal
#' @noRd
mc_trend_reference <- function(levs, reference = NULL) {
    if (is.null(reference)) {
        # Two groups is nearly always treatment against control, so the
        # deviation is reported from a level. Beyond two there is no
        # defensible alphabetical baseline, so every group is reported
        # against the trend over all of them instead.
        if (length(levs) == 2L) {
            return(list(
                group = "ordered", reference = levs[1], levels = levs,
                set_by = "alphabetical default"
            ))
        }
        return(list(
            group = "global", reference = "global", levels = levs,
            set_by = "default for more than two groups"
        ))
    }

    if (length(reference) != 1L || !is.character(reference)) {
        stop(
            "reference must be a single level name, or \"global\".",
            call. = FALSE
        )
    }
    if (identical(reference, "global")) {
        return(list(
            group = "global", reference = "global", levels = levs,
            set_by = "set by user"
        ))
    }
    if (!reference %in% levs) {
        stop(
            "reference must be one level of the group column, or \"global\". ",
            "You gave '", reference, "'; available: ",
            paste(levs, collapse = ", "), ".",
            call. = FALSE
        )
    }
    list(
        group = "ordered", reference = reference,
        levels = c(reference, setdiff(levs, reference)),
        set_by = "set by user"
    )
}


#' One empty result row
#'
#' @description
#' A taxon that could not be fitted still gets a row, so that the output has
#' one row per taxon asked for and a reader can see which failed and why.
#'
#' @param taxon The taxon name.
#' @param n_obs,n_subjects What the taxon had.
#' @param note Why it failed, or `NA`.
#'
#' @return A one-row data.frame.
#'
#' @keywords internal
#' @noRd
mc_trend_blank <- function(taxon, n_obs = 0L, n_subjects = 0L, note = NA_character_) {
    data.frame(
        taxon = taxon, edf = NA_real_, F = NA_real_, p_value = NA_real_,
        n_obs = as.integer(n_obs), n_subjects = as.integer(n_subjects),
        converged = FALSE, note = note, stringsAsFactors = FALSE
    )
}


#' Fit one taxon and read the test off the fit
#'
#' @param d The taxon's data: `value`, `time`, `subject`, covariates, and
#'   `mc_grp` when there is a group.
#' @param k Basis size.
#' @param subject_effect The repeated-measures term.
#' @param covariates Character vector, or `NULL`.
#' @param group One of `"none"`, `"ordered"`, `"global"`.
#' @param taxon The taxon's name, for the returned row.
#'
#' @return A one-row data.frame of `taxon`, `edf`, `F`, `p_value`, `n_obs`,
#'   `n_subjects`, `converged`, `note`.
#'
#' @keywords internal
#' @noRd
mc_trend_fit_one <- function(d, k, subject_effect, covariates, group, taxon) {
    d <- d[!is.na(d$value), , drop = FALSE]
    n_obs <- nrow(d)
    n_sub <- length(unique(d$subject))

    # A taxon seen at fewer distinct times than the basis has functions
    # cannot support the smooth, and a group that lost all of its samples
    # leaves a factor level with no data.
    n_t <- length(unique(d$time))
    if (n_obs < k + 2L || n_t < 3L) {
        return(mc_trend_blank(taxon, n_obs, n_sub, "too few observations"))
    }
    if (!identical(group, "none")) {
        d$mc_grp <- droplevels(d$mc_grp)
        if (nlevels(d$mc_grp) < 2L) {
            return(mc_trend_blank(taxon, n_obs, n_sub, "fewer than 2 groups"))
        }
    }
    if (stats::var(d$value) < .Machine$double.eps) {
        return(mc_trend_blank(taxon, n_obs, n_sub, "constant"))
    }

    d$subject <- factor(d$subject)
    kk <- min(k, n_t - 1L)

    out <- tryCatch(
        {
            f1 <- mc_trend_formula(
                mc_trend_terms(kk, subject_effect, covariates, group, TRUE)
            )
            m1 <- mc_trend_gam(f1, d)

            if (identical(group, "none")) {
                mc_trend_read_smooth(m1, "s(time)")
            } else if (nlevels(d$mc_grp) == 2L) {
                # Two groups leave exactly one difference smooth, so its own
                # p-value is already the whole test and the null model need
                # not be fitted. The test is chosen by the number of groups
                # and not by the parameterisation: with a named reference and
                # three groups there are several difference smooths, and
                # reporting whichever came first would be the arbitrary
                # choice the omnibus exists to avoid.
                mc_trend_read_smooth(m1, NULL)
            } else {
                f0 <- mc_trend_formula(
                    mc_trend_terms(kk, subject_effect, covariates, group, FALSE)
                )
                m0 <- mc_trend_gam(f0, d)
                mc_trend_read_anova(m0, m1)
            }
        },
        error = function(e) {
            list(
                edf = NA_real_, F = NA_real_, p = NA_real_, ok = FALSE,
                note = sub("[\r\n]+", " ", conditionMessage(e))
            )
        }
    )

    data.frame(
        taxon = taxon, edf = out$edf, F = out$F, p_value = out$p,
        n_obs = as.integer(n_obs), n_subjects = as.integer(n_sub),
        converged = isTRUE(out$ok), note = out$note %||% NA_character_,
        stringsAsFactors = FALSE
    )
}


#' Fit a gam, muffling one expected diagnostic
#'
#' @description
#' `subject_effect = "slope"` and `"smooth"` both put a second smooth on a
#' variable that `s(time)` or the subject term already uses, and `mgcv`
#' reports that as repeated one-dimensional smooths of the same variable.
#' For these structures it is expected rather than a problem: the factor
#' smooth is built with `m = 1`, which penalises its null space so that the
#' global smooth stays identifiable, and the random slope is a product basis
#' rather than a competing smooth. Left alone the message would be raised
#' once per taxon, so this one message is muffled and every other warning is
#' passed through untouched.
#'
#' @param f The formula.
#' @param d The taxon's data.
#'
#' @return A fitted `gam`.
#'
#' @keywords internal
#' @noRd
mc_trend_gam <- function(f, d) {
    withCallingHandlers(
        mgcv::gam(f, data = d, method = "REML"),
        warning = function(w) {
            if (grepl("repeated 1-d smooths", conditionMessage(w))) {
                invokeRestart("muffleWarning")
            }
        }
    )
}


#' Read a single smooth's test out of summary.gam
#'
#' @param m A fitted `gam`.
#' @param want The smooth's label, or `NULL` to take the `by =` smooths.
#'
#' @return A list of `edf`, `F`, `p`, `ok` and `note`.
#'
#' @keywords internal
#' @noRd
mc_trend_read_smooth <- function(m, want = NULL) {
    st <- summary(m)$s.table
    if (is.null(st) || nrow(st) == 0) {
        return(list(
            edf = NA_real_, F = NA_real_, p = NA_real_, ok = FALSE,
            note = "no smooth terms in the fit"
        ))
    }
    rn <- rownames(st)

    # The difference smooth is named for the level it belongs to, which
    # varies, so it is found by the colon that only a by-smooth carries.
    rows <- if (is.null(want)) {
        which(grepl("):", rn, fixed = TRUE))
    } else {
        which(rn == want)
    }
    if (length(rows) == 0) {
        return(list(
            edf = NA_real_, F = NA_real_, p = NA_real_, ok = FALSE,
            note = "expected smooth term absent from the fit"
        ))
    }
    i <- rows[1]
    pcol <- grep("^p-value$|^Pr", colnames(st))[1]
    fcol <- grep("^F$|^Chi", colnames(st))[1]
    list(
        edf = unname(st[i, "edf"]),
        F = if (is.na(fcol)) NA_real_ else unname(st[i, fcol]),
        p = if (is.na(pcol)) NA_real_ else unname(st[i, pcol]),
        ok = isTRUE(m$converged) || is.null(m$converged),
        note = NA_character_
    )
}


#' Read the nested comparison of two fits
#'
#' @description
#' The null model holds one curve for every group, the alternative lets each
#' group depart from it. The comparison is therefore the omnibus question
#' "does trajectory shape depend on group", with one answer however many
#' groups there are.
#'
#' The p-value is approximate: the smoothing parameters were estimated rather
#' than known, which `mgcv::anova.gam` documents as a limitation of the test.
#'
#' One consequence of that has to be handled here. Each model chooses its own
#' smoothing parameters, so adding the group smooths can leave the deviance
#' no lower than the null model's. `mgcv` then cannot form the F statistic
#' and returns nothing. That is not a failed fit: it says the group smooths
#' bought no improvement, which is the strongest possible absence of
#' evidence that shape depends on group, so it is reported as `p = 1` with a
#' note rather than as a missing result that would drop the taxon out of the
#' correction across taxa.
#'
#' @param m0,m1 Fitted `gam` objects, `m0` nested in `m1`.
#'
#' @return A list of `edf`, `F`, `p`, `ok` and `note`.
#'
#' @keywords internal
#' @noRd
mc_trend_read_anova <- function(m0, m1) {
    a <- suppressWarnings(stats::anova(m0, m1, test = "F"))
    pcol <- grep("^Pr", colnames(a))[1]
    fcol <- grep("^F$|^Chi", colnames(a))[1]

    pv <- if (is.na(pcol)) NA_real_ else a[[pcol]]
    fv <- if (is.na(fcol)) NA_real_ else a[[fcol]]
    pv <- pv[!is.na(pv)]
    fv <- fv[!is.na(fv)]

    st <- summary(m1)$s.table
    by_rows <- if (is.null(st)) integer(0) else {
        which(grepl("):", rownames(st), fixed = TRUE))
    }
    edf <- if (length(by_rows) > 0) sum(st[by_rows, "edf"]) else NA_real_
    ok <- (isTRUE(m1$converged) || is.null(m1$converged)) &&
        (isTRUE(m0$converged) || is.null(m0$converged))

    if (length(pv) == 0 && isTRUE(ok)) {
        dev <- a[["Deviance"]]
        dev <- dev[!is.na(dev)]
        # The column is the null model's residual deviance minus the
        # alternative's, so a value at or below zero means the extra terms
        # did not improve the fit.
        if (length(dev) > 0 && dev[length(dev)] <= 0) {
            return(list(
                edf = edf, F = NA_real_, p = 1, ok = TRUE,
                note = "group smooths did not reduce the deviance"
            ))
        }
    }

    list(
        edf = edf,
        F = if (length(fv) > 0) fv[length(fv)] else NA_real_,
        p = if (length(pv) > 0) pv[length(pv)] else NA_real_,
        ok = ok,
        note = if (length(pv) > 0) {
            NA_character_
        } else {
            "the nested comparison yielded no F statistic"
        }
    )
}


#' Run a per-taxon function over every taxon
#'
#' @description
#' Taxa are independent of one another, so this is the level the work is
#' split at. `parallel::mclapply` is used where the platform forks; on
#' Windows it does not, and the run is serial with one message saying so
#' rather than a cluster whose setup cost would exceed the saving on the
#' taxon counts this is used at.
#'
#' @param n Number of taxa.
#' @param fun Function of one index.
#' @param n_cores Integer.
#' @param verbose Logical.
#'
#' @return A list of results.
#'
#' @keywords internal
#' @noRd
mc_trend_apply <- function(n, fun, n_cores = 1L, verbose = TRUE) {
    n_cores <- as.integer(n_cores)
    if (is.na(n_cores) || n_cores < 1L) n_cores <- 1L

    if (n_cores > 1L && .Platform$OS.type != "unix") {
        if (isTRUE(verbose)) {
            message(
                "  n_cores = ", n_cores, " ignored: forking is not available ",
                "on this platform, running on one core"
            )
        }
        n_cores <- 1L
    }
    if (n_cores <= 1L) {
        return(lapply(seq_len(n), fun))
    }
    parallel::mclapply(
        seq_len(n), fun,
        mc.cores = min(n_cores, n), mc.preschedule = FALSE
    )
}


#' Adjust p-values across taxa, leaving the failures as NA
#'
#' @param p Numeric vector of p-values.
#' @param method A `stats::p.adjust` method.
#'
#' @return A numeric vector the same length as `p`.
#'
#' @keywords internal
#' @noRd
mc_trend_adjust <- function(p, method = "BH") {
    if (!method %in% stats::p.adjust.methods) {
        stop(
            "p_adjust must be one of ",
            mc_fmt_some(stats::p.adjust.methods), ".",
            call. = FALSE
        )
    }
    q <- rep(NA_real_, length(p))
    ok <- !is.na(p)
    # A taxon that failed to fit is not a test that was carried out, so it
    # is left out of the correction rather than counted in its denominator.
    if (any(ok)) q[ok] <- stats::p.adjust(p[ok], method = method)
    q
}


#' Warn when a per-subject shape can swallow the group effect
#'
#' @description
#' If every sample from a subject carries the same group, then letting each
#' subject have its own curve lets the subject curves absorb the difference
#' between groups. The fit still runs; the test loses power. The caller is
#' told once rather than left to wonder why nothing reached significance.
#'
#' @param design The aligned design.
#' @param subject_effect The chosen term.
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
#' @noRd
mc_trend_warn_absorbed <- function(design, subject_effect) {
    if (!identical(subject_effect, "smooth") || is.null(design$mc_grp)) {
        return(invisible(NULL))
    }
    per <- tapply(design$mc_grp, design$subject, function(g) {
        length(unique(g))
    })
    if (all(per == 1L)) {
        warning(
            "subject_effect = \"smooth\" with a group that is constant ",
            "within each subject: the per-subject curves can absorb the ",
            "group difference being tested, which lowers power. Compare ",
            "against subject_effect = \"intercept\".",
            call. = FALSE
        )
    }
    invisible(NULL)
}


#' Write a trend or comparison result to disk
#'
#' @param res The results data.frame.
#' @param meta Named list of lines describing the run.
#' @param out_dir Directory to write into.
#' @param stem File name stem, without extension.
#'
#' @return Character vector of the paths written.
#'
#' @keywords internal
#' @noRd
mc_trend_write <- function(res, meta, out_dir, stem) {
    if (!dir.exists(out_dir)) {
        ok <- dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        if (!ok) {
            stop("Could not create out_dir '", out_dir, "'.", call. = FALSE)
        }
    }
    f <- file.path(out_dir, paste0(stem, ".tsv"))
    utils::write.table(
        res, f,
        sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
    )
    log <- file.path(out_dir, paste0(stem, "_log.txt"))
    writeLines(meta, log)
    c(f, log)
}
