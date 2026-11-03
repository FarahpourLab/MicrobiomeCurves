############################################
# Trend testing.
#
# Two questions, two functions:
#
#   mc_trends()   does this taxon change over time at all?
#   mc_compare()  do groups differ in the shape of their trajectory?
#
# Both fit a generalised additive model per taxon and differ only in the
# group term, so the shared machinery lives in trends_helpers.R. Neither
# touches the imputation path: a completed table from mc_run() is one valid
# input and a table with gaps is another, because a GAM drops incomplete
# rows and fits on what is left.
############################################


#' Test whether a taxon changes over time
#'
#' @description
#' Fits, for every taxon, a smooth of abundance against time with a
#' penalised random effect for subject, and reports whether the trajectory
#' departs from flat.
#'
#' @details
#' The model is
#'
#' \preformatted{
#'   value ~ s(time, k = k) + <covariates> + <subject term>
#' }
#'
#' fitted by `mgcv::gam()` with `method = "REML"`. The reported test is the
#' p-value of `s(time)` from `summary.gam()`, adjusted across taxa.
#'
#' Time must be numeric. A smooth needs a continuous axis, so a column of
#' labels such as `"week_4"` is refused rather than coerced to 1, 2, 3,
#' which would assert an equal spacing the study may not have.
#'
#' `k` is the size of the spline basis: how much shape the curve is allowed
#' to have. `mgcv` defaults to 10, which fails below ten distinct time
#' points, so the default here is `min(n_distinct_times - 1, 10)`, read from
#' the design. The penalty, not `k`, decides how wiggly the fitted curve
#' actually is, which is what the reported `edf` measures: `edf` near 1 means
#' the fit was flattened to a line.
#'
#' `subject_effect` says how repeated measures from one subject are handled.
#' Treating them as independent would be pseudo-replication and would make
#' every taxon look significant.
#'
#' \itemize{
#'   \item `"intercept"` — `s(subject, bs = "re")`. Subjects differ in
#'     level but share a shape. The usual case in microbiome data, and the
#'     default.
#'   \item `"slope"` — adds `s(subject, time, bs = "re")`. Subjects also
#'     drift apart over time, so between-subject spread grows. Needs about
#'     four time points per subject.
#'   \item `"smooth"` — `s(time, subject, bs = "fs", m = 1)`. Each subject
#'     gets its own curve, penalised toward the mean. Needs about six time
#'     points per subject, and in [mc_compare()] it can absorb the group
#'     difference being tested; a warning is raised when that applies.
#'   \item `"none"` — nothing added. For one sample per subject.
#' }
#'
#' These structures are the ones described for `mgcv` by Wood (2017) and, for
#' group-varying smooths, by Pedersen et al. (2019).
#'
#' @param abundance Abundance table: a data.frame with one column of taxon
#'   names and one column per sample, or a matrix of taxa by samples with
#'   row names. Values are expected on the centred log ratio scale. `NA` is
#'   allowed and those rows are dropped from the fit, so passing a table
#'   with gaps tests observed cells only while passing `mc_run()`'s
#'   `completed` table includes the imputed ones.
#' @param metadata A data.frame with one row per sample.
#' @param sample_col Column of `metadata` holding the sample names used as
#'   `abundance`'s column names.
#' @param subject_col Column of `metadata` saying which subject a sample
#'   came from.
#' @param time_col Column of `metadata` holding numeric time.
#' @param taxon_col Column of `abundance` holding the taxon names. `NULL`
#'   takes the first column that does not hold numbers.
#' @param normalization What to do to each sample before fitting. `"NONE"`
#'   leaves the table alone and is the default, because `mc_run()`'s
#'   `completed` table is already on the centred log ratio scale. `"TSS"`
#'   divides each sample by its total. `"CLR"` takes centred log ratios by
#'   the same step `mc_run(abundance_type = "raw")` applies, so an
#'   untransformed abundance table can be tested directly and a trend test
#'   agrees value for value with an imputation of that table. Applied to the whole table before
#'   `taxa` selects from it, since both normalisations are properties of
#'   the sample as a whole. The names follow MaAsLin 2.
#' @param transform What to do to each value after normalisation. `"NONE"`
#'   is the default and leaves the values alone. The others are `"LOG"`
#'   (log2, with zeros replaced by half the smallest positive value the
#'   taxon reaches), `"LOGIT"`, and `"AST"`, the arcsine square root.
#'   `"CLR"` normalisation leaves values on the whole real line, so it
#'   accepts only `"NONE"`; any other pairing stops the call.
#' @param pseudocount Used only by `normalization = "CLR"`, and meaning
#'   exactly what it means in [mc_run()]: `"auto"` replaces each zero with
#'   0.65 of the smallest positive value in its own sample, and a number is
#'   added to every entry instead. `transform = "LOG"` does not read it, so
#'   that `"TSS"` with `"LOG"` reproduces MaAsLin 2's own rule of half the
#'   smallest positive value in the taxon.
#' @param covariates Character vector of `metadata` columns to adjust for.
#'   They enter as parametric terms, shifting the curve without a shape of
#'   their own, and no test is reported for them.
#' @param k Spline basis size. `NULL`, the default, takes the largest the
#'   design supports: one fewer than the number of distinct time points. A
#'   single number sets it. Several numbers are candidates rather than a
#'   size, and each taxon keeps the one `k_select` prefers, capped by the
#'   time points that taxon actually has. The chosen size is reported in
#'   the `k` column of `results`.
#' @param k_select How a candidate `k` is chosen, used only when `k` holds
#'   more than one value. `"aic"` takes the lowest AIC. `"loo"` takes the
#'   lowest leave-one-out mean squared error, computed in closed form from
#'   the leverages rather than by refitting; note that it leaves out one
#'   sample and not one subject, so with repeated measures it flatters the
#'   larger bases. `"variance"` takes the smallest candidate whose deviance
#'   explained reaches `k_var` of the largest candidate's, which asks for
#'   the smallest basis that loses nothing worth keeping rather than for
#'   the best-fitting one.
#' @param k_var The fraction used by `k_select = "variance"`, between 0 and
#'   1. It is a fraction of what the largest candidate explains, not an
#'   absolute deviance explained: few microbiome taxa reach 0.99 of their
#'   own variance, and a threshold none of them clears would choose nothing.
#' @param subject_effect One of `"intercept"`, `"slope"`, `"smooth"` or
#'   `"none"`.
#' @param taxa Character vector of taxa to test, or `NULL` for all.
#' @param p_adjust A `stats::p.adjust` method.
#' @param alpha Significance threshold on the adjusted p-value. It is
#'   reported, and [mc_trend_plot()] takes it as the default set of taxa to
#'   draw; it does not affect any fit.
#' @param out_dir Directory to write the results table into, or `NULL` to
#'   write nothing. No figures are written: drawing is [mc_trend_plot()]'s
#'   job, so the table can be looked over before anything is drawn.
#' @param n_cores Number of cores. Taxa are independent, so this splits by
#'   taxon. Forking only; ignored on Windows.
#' @param verbose Logical. Whether to report the resolved settings.
#'
#' @return An object of class `mc_trends` with
#' \itemize{
#'   \item `results`: one row per taxon, ordered by `q_value`. The `k`
#'     column holds the basis size that taxon was fitted with.
#'   \item `design`, `abundance`: what was fitted, kept so a follow-up can
#'     refit a taxon without the models being stored.
#'   \item `settings`: the resolved `k`, subject term and covariates.
#'   \item `files`: anything written.
#' }
#'
#' @references
#' Wood SN (2017). *Generalized Additive Models: An Introduction with R*,
#' 2nd edition. Chapman and Hall/CRC.
#'
#' Pedersen EJ, Miller DL, Simpson GL, Ross N (2019). Hierarchical
#' generalized additive models in ecology: an introduction with mgcv.
#' *PeerJ* 7:e6876.
#'
#' @examples
#' set.seed(1)
#' meta <- expand.grid(
#'     time = 1:6, subject = paste0("S", 1:12),
#'     stringsAsFactors = FALSE
#' )
#' meta$sample <- paste0(meta$subject, "_T", meta$time)
#'
#' # one taxon with a trend, one without
#' m <- matrix(rnorm(2 * nrow(meta)), nrow = 2,
#'             dimnames = list(c("Rising", "Flat"), meta$sample))
#' m["Rising", ] <- m["Rising", ] + 0.8 * meta$time
#' tab <- data.frame(taxon = rownames(m), m, check.names = FALSE)
#'
#' tr <- mc_trends(tab, meta, "sample", "subject", "time", verbose = FALSE)
#' tr$results[, c("taxon", "edf", "p_value", "q_value")]
#'
#' @seealso [mc_trend_plot()] to draw the result, [mc_compare()] to ask
#'   whether groups differ, [mc_run()] to impute missing samples first.
#'
#' @importFrom mgcv gam s
#' @export
mc_trends <- function(abundance, metadata,
                      sample_col, subject_col, time_col,
                      taxon_col = NULL,
                      normalization = c("NONE", "TSS", "CLR"),
                      transform = c("NONE", "LOG", "LOGIT", "AST"),
                      pseudocount = "auto",
                      covariates = NULL,
                      k = NULL,
                      k_select = c("aic", "loo", "variance"),
                      k_var = 0.99,
                      subject_effect = c("intercept", "slope", "smooth", "none"),
                      taxa = NULL,
                      p_adjust = "BH",
                      alpha = 0.05,
                      out_dir = NULL,
                      n_cores = 1,
                      verbose = TRUE) {
    subject_effect <- match.arg(subject_effect)
    normalization <- match.arg(normalization)
    transform <- match.arg(transform)
    k_select <- match.arg(k_select)
    mc_trend_check_k_var(k_var)

    mc_trend_run(
        abundance = abundance, metadata = metadata,
        sample_col = sample_col, subject_col = subject_col,
        time_col = time_col, taxon_col = taxon_col,
        normalization = normalization, transform = transform,
        pseudocount = pseudocount,
        group_col = NULL, reference = NULL,
        covariates = covariates, k = k, k_select = k_select,
        k_var = k_var, subject_effect = subject_effect,
        taxa = taxa, p_adjust = p_adjust, alpha = alpha,
        out_dir = out_dir, n_cores = n_cores, verbose = verbose
    )
}


#' Test whether groups differ in the shape of their trajectory
#'
#' @description
#' Fits, for every taxon, one model holding both the overall smooth of time
#' and a smooth per group, and reports whether the groups' trajectories
#' differ in shape.
#'
#' @details
#' One model is fitted, not one per group. A separate fit per group would
#' give two curves but no standard error for the gap between them, and could
#' never produce "the groups are the same", because two independently fitted
#' curves always differ somewhere. In the joint model the difference is an
#' estimated term with a penalty of its own, so "same shape" is the penalty
#' shrinking that term to zero, and that is what gets tested.
#'
#' This is the opposite of what [mc_run()] does with `group_col`, and
#' deliberately: imputation must not borrow a trajectory across groups, or
#' the group effect leaks into the imputed values, while comparing requires
#' both groups in one fit.
#'
#' Two parameterisations, chosen by `reference`:
#'
#' \itemize{
#'   \item A level name gives an ordered factor, so `s(time)` is that
#'     level's own curve and the `by =` smooth is the deviation from it. With
#'     two groups that reads as "treated minus control".
#'   \item `"global"` gives an unordered factor, so the `by =` smooths are
#'     centred, `s(time)` becomes the trend over all groups, and each smooth
#'     is one group's departure from it. No baseline has to be chosen and
#'     every group gets a departure of its own.
#' }
#'
#' The parameterisation decides what the fitted terms mean; the **number of
#' groups** decides how the single reported p-value is obtained. Two groups
#' leave exactly one difference smooth, whose own p-value from
#' `summary.gam()` is already the whole test, so one model is fitted. Three
#' or more leave several, and reporting one of them would depend on which
#' level happened to come first, so the null model is fitted too and the
#' answer is the omnibus `anova(m0, m1, test = "F")`: one p-value meaning
#' "does trajectory shape depend on group", for any number of groups and
#' whichever parameterisation is in use.
#'
#' `reference = NULL` takes the first level alphabetically when there are two
#' groups, since that is almost always treatment against control, and
#' `"global"` beyond two, since no alphabetical baseline among several arms
#' is defensible. Whichever is used is reported when `verbose`.
#'
#' Note on the nested test: `mgcv` documents the p-value from comparing two
#' fits whose smoothing parameters were estimated rather than known as
#' approximate. It is used here because the models differ by whole terms,
#' but it is not exact.
#'
#' Each group's `s(time)` is estimated with all the data contributing, not
#' from that group's samples alone, so a reported reference curve is
#' regularised toward the pooled one.
#'
#' @inheritParams mc_trends
#' @param group_col Column of `metadata` holding the variable being tested.
#'   Required. Its levels are allowed their own shape over time, and the
#'   reported p-value is about it.
#' @param reference A level of `group_col` to measure deviations from, or
#'   `"global"` to measure them from the trend over all groups. `NULL`
#'   resolves by the number of levels.
#'
#' @return An object of class `mc_compare`, as [mc_trends()] returns, with
#'   `n_groups` and `reference` added to `results` and the resolved group
#'   settings in `settings`.
#'
#' @references
#' Wood SN (2017). *Generalized Additive Models: An Introduction with R*,
#' 2nd edition. Chapman and Hall/CRC.
#'
#' Pedersen EJ, Miller DL, Simpson GL, Ross N (2019). Hierarchical
#' generalized additive models in ecology: an introduction with mgcv.
#' *PeerJ* 7:e6876.
#'
#' @examples
#' set.seed(1)
#' meta <- expand.grid(
#'     time = 1:6, subject = paste0("S", 1:12),
#'     stringsAsFactors = FALSE
#' )
#' meta$sample <- paste0(meta$subject, "_T", meta$time)
#' meta$arm <- ifelse(
#'     meta$subject %in% paste0("S", 1:6), "control", "treated"
#' )
#'
#' # one taxon where the arms diverge, one where they do not
#' m <- matrix(rnorm(2 * nrow(meta)), nrow = 2,
#'             dimnames = list(c("Diverges", "Same"), meta$sample))
#' m["Diverges", ] <- m["Diverges", ] +
#'     ifelse(meta$arm == "treated", 0.9 * meta$time, 0)
#' tab <- data.frame(taxon = rownames(m), m, check.names = FALSE)
#'
#' cm <- mc_compare(
#'     tab, meta, "sample", "subject", "time",
#'     group_col = "arm", verbose = FALSE
#' )
#' cm$results[, c("taxon", "edf", "p_value", "q_value", "reference")]
#'
#' @seealso [mc_trends()] to ask only about time.
#'
#' @importFrom mgcv gam s
#' @export
mc_compare <- function(abundance, metadata,
                       sample_col, subject_col, time_col,
                       group_col,
                       reference = NULL,
                       taxon_col = NULL,
                       normalization = c("NONE", "TSS", "CLR"),
                       transform = c("NONE", "LOG", "LOGIT", "AST"),
                       pseudocount = "auto",
                       covariates = NULL,
                       k = NULL,
                       k_select = c("aic", "loo", "variance"),
                       k_var = 0.99,
                       subject_effect = c("intercept", "slope", "smooth", "none"),
                       taxa = NULL,
                       p_adjust = "BH",
                       alpha = 0.05,
                       out_dir = NULL,
                       n_cores = 1,
                       verbose = TRUE) {
    subject_effect <- match.arg(subject_effect)
    normalization <- match.arg(normalization)
    transform <- match.arg(transform)
    k_select <- match.arg(k_select)
    mc_trend_check_k_var(k_var)

    if (missing(group_col) || is.null(group_col)) {
        stop(
            "group_col is required: mc_compare() asks whether groups differ. ",
            "Use mc_trends() to ask only whether a taxon changes over time.",
            call. = FALSE
        )
    }
    if (length(group_col) != 1L || !is.character(group_col)) {
        stop("group_col must be a single column name.", call. = FALSE)
    }

    mc_trend_run(
        abundance = abundance, metadata = metadata,
        sample_col = sample_col, subject_col = subject_col,
        time_col = time_col, taxon_col = taxon_col,
        normalization = normalization, transform = transform,
        pseudocount = pseudocount,
        group_col = group_col, reference = reference,
        covariates = covariates, k = k, k_select = k_select,
        k_var = k_var, subject_effect = subject_effect,
        taxa = taxa, p_adjust = p_adjust, alpha = alpha,
        out_dir = out_dir, n_cores = n_cores, verbose = verbose
    )
}


#' The shared body of mc_trends() and mc_compare()
#'
#' @description
#' The two exported functions differ in one term of the model and in two
#' arguments. Everything else - validation, alignment, the per-taxon loop,
#' the correction across taxa, the output - is the same, and is here.
#'
#' @inheritParams mc_compare
#'
#' @return An `mc_trends` or `mc_compare` object.
#'
#' @keywords internal
#' @noRd
mc_trend_run <- function(abundance, metadata, sample_col, subject_col,
                         time_col, taxon_col, normalization, transform,
                         pseudocount, group_col, reference,
                         covariates, k, k_select, k_var, subject_effect,
                         taxa, p_adjust,
                         alpha, out_dir, n_cores, verbose) {
    say <- function(...) if (isTRUE(verbose)) message(...)

    design <- mc_trend_design(
        metadata, sample_col, subject_col, time_col, group_col, covariates
    )
    tm <- mc_trend_matrix(
        abundance, taxon_col, taxa,
        normalization, transform, pseudocount
    )
    al <- mc_trend_align(design, tm$mat, verbose)
    design <- al$design
    mat <- al$mat

    k <- mc_trend_k(design$time, k)

    grp <- list(group = "none", reference = NA_character_, levels = NULL)
    if (!is.null(group_col)) {
        levs <- sort(unique(design$mc_grp))
        if (length(levs) < 2L) {
            stop(
                "group_col '", group_col, "' has ", length(levs),
                " level(s) in the samples being used. At least 2 are needed.",
                call. = FALSE
            )
        }
        grp <- mc_trend_reference(levs, reference)
        design <- mc_trend_group_factor(design, grp$levels, grp$group)
    }
    mc_trend_warn_absorbed(design, subject_effect)

    n_tax <- nrow(mat)
    say(
        if (is.null(group_col)) "mc_trends: " else "mc_compare: ",
        n_tax, " taxa, ", nrow(design), " samples, ",
        length(unique(design$subject)), " subjects, ",
        length(unique(design$time)), " time points"
    )
    say("  k              : ", if (length(k) == 1L) {
        k
    } else {
        paste0(
            mc_fmt_some(k), ", chosen per taxon by ",
            if (identical(k_select, "aic")) {
                "AIC"
            } else if (identical(k_select, "loo")) {
                "leave-one-out MSE"
            } else {
                paste0(k_var, " of the largest candidate's deviance")
            }
        )
    })
    say("  subject effect : ", subject_effect)
    # Reported only when something was applied. Silence is the default case
    # and saying so on every run would bury the lines that vary.
    if (!identical(normalization, "NONE") || !identical(transform, "NONE")) {
        for (ln in mc_trend_scale_lines(
            normalization, transform, pseudocount
        )) {
            say(ln)
        }
    }
    if (length(covariates) > 0) {
        say("  covariates     : ", mc_fmt_some(covariates))
    }
    if (!is.null(group_col)) {
        say("  group          : ", group_col, " (", length(grp$levels),
            " levels: ", mc_fmt_some(grp$levels, n = 8), ")")
        say("  reference      : ", grp$reference, "  [", grp$set_by, "]")
        say("  reporting      : ", if (identical(grp$group, "global")) {
            "each group minus the trend over all groups"
        } else {
            paste0("each group minus ", grp$reference)
        })
        say("  group sizes    : ", paste(
            names(table(design$mc_grp)), as.integer(table(design$mc_grp)),
            sep = " ", collapse = ", "
        ))
    }

    one <- function(i) {
        d <- design
        d$value <- as.numeric(mat[i, ])
        mc_trend_fit_one(
            d, k, subject_effect, covariates, grp$group, rownames(mat)[i],
            k_select, k_var
        )
    }
    res <- do.call(rbind, mc_trend_apply(n_tax, one, n_cores, verbose))

    res$q_value <- mc_trend_adjust(res$p_value, p_adjust)
    if (!is.null(group_col)) {
        res$n_groups <- length(grp$levels)
        res$reference <- grp$reference
    }
    front <- c("taxon", "edf", "F", "p_value", "q_value")
    res <- res[, c(front, setdiff(names(res), front)), drop = FALSE]
    res <- res[order(res$q_value, res$p_value, res$taxon,
        na.last = TRUE
    ), , drop = FALSE]
    rownames(res) <- NULL

    n_fail <- sum(!res$converged)
    n_hit <- sum(!is.na(res$q_value) & res$q_value < alpha)
    say("  tested         : ", n_tax - n_fail, " of ", n_tax,
        if (n_fail > 0) paste0(" (", n_fail, " failed)") else "")
    say("  q < ", alpha, "        : ", n_hit)

    obj <- structure(
        list(
            results = res,
            design = design,
            abundance = mat,
            settings = list(
                k = k, k_select = k_select, k_var = k_var,
                subject_effect = subject_effect,
                covariates = covariates, group_col = group_col,
                normalization = normalization, transform = transform,
                pseudocount = pseudocount,
                reference = grp$reference, group = grp$group,
                levels = grp$levels, alpha = alpha, p_adjust = p_adjust,
                sample_col = sample_col, subject_col = subject_col,
                time_col = time_col,
                # What was fitted and what was tested, kept here because a
                # result that cannot say which model produced it cannot be
                # reproduced, and because nothing else records the choice:
                # the number of group levels decides the test, and that is
                # not visible from the table of p-values.
                # With candidates there is no one k to print, and printing
                # any single one of them would misreport every taxon that
                # chose another. The per-taxon sizes are in results$k.
                formula = mc_trend_formula_text(
                    k, subject_effect, covariates, grp$group
                ),
                method = "REML",
                test = if (identical(grp$group, "none")) {
                    "s(time) from summary.gam()"
                } else if (length(grp$levels) == 2L) {
                    "difference smooth from summary.gam()"
                } else {
                    "anova(m0, m1, test = \"F\"), omnibus over groups"
                }
            ),
            files = character(0)
        ),
        class = if (is.null(group_col)) "mc_trends" else "mc_compare"
    )

    files <- character(0)
    stem <- if (is.null(group_col)) "trends" else "compare"
    if (!is.null(out_dir)) {
        files <- mc_trend_write(res, out_dir, stem)
    }
    obj$files <- files
    if (length(files) > 0) {
        say("  wrote ", length(files), " file(s) to ", out_dir)
    }
    obj
}


#' The lines of a trend run log
#'
#' @param obj The result object.
#' @param n_tax,n_fail,n_hit Counts from the run.
#'
#' @return A character vector of lines.
#'
#' @keywords internal
#' @noRd
mc_trend_log_lines <- function(obj, n_tax, n_fail, n_hit) {
    s <- obj$settings
    d <- obj$design
    grouped <- !is.null(s$group_col)

    terms <- mc_trend_terms(
        s$k, s$subject_effect, s$covariates, s$group, TRUE
    )

    c(
        paste0(
            "MicrobiomeCurves ",
            if (grouped) "mc_compare" else "mc_trends", " log"
        ),
        paste0("package version: ", mc_version_string()),
        "",
        "DESIGN",
        paste0("  taxa        : ", n_tax),
        paste0("  samples     : ", nrow(d)),
        paste0("  subjects    : ", length(unique(d$subject))),
        paste0("  time points : ", length(unique(d$time))),
        paste0("  time range  : ", min(d$time), " to ", max(d$time)),
        "",
        "INPUT SCALE",
        mc_trend_scale_lines(s$normalization, s$transform, s$pseudocount),
        "",
        "MODEL",
        paste0("  value ~ ", paste(terms, collapse = " + ")),
        paste0("  fitted by mgcv::gam(method = \"REML\")"),
        paste0("  k              : ", s$k),
        paste0("  subject effect : ", s$subject_effect),
        paste0(
            "  covariates     : ",
            if (length(s$covariates) > 0) mc_fmt_some(s$covariates) else "none"
        ),
        if (grouped) {
            c(
                paste0("  group column   : ", s$group_col),
                paste0("  group levels   : ", mc_fmt_some(s$levels, n = 20)),
                paste0("  reference      : ", s$reference),
                paste0(
                    "  test           : ",
                    if (length(s$levels) == 2L) {
                        "difference smooth from summary.gam()"
                    } else {
                        "anova(m0, m1, test = \"F\"), omnibus over groups"
                    }
                ),
                paste0(
                    "  note           : a p-value from comparing fits whose ",
                    "smoothing"
                ),
                "                   parameters were estimated is approximate."
            )
        } else {
            paste0("  test           : s(time) from summary.gam()")
        },
        "",
        "RESULT",
        paste0("  tested      : ", n_tax - n_fail, " of ", n_tax),
        paste0("  failed      : ", n_fail),
        paste0("  p_adjust    : ", s$p_adjust),
        paste0("  q < ", s$alpha, "   : ", n_hit)
    )
}


#' Print a time-trend result
#'
#' @param x An `mc_trends` object, as returned by [mc_trends()].
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#'
#' @examples
#' set.seed(1)
#' meta <- expand.grid(
#'     time = 1:6, subject = paste0("S", 1:8),
#'     stringsAsFactors = FALSE
#' )
#' meta$sample <- paste0(meta$subject, "_T", meta$time)
#' m <- matrix(rnorm(nrow(meta)), nrow = 1,
#'             dimnames = list("Rising", meta$sample))
#' m[1, ] <- m[1, ] + 0.8 * meta$time
#' tab <- data.frame(taxon = rownames(m), m, check.names = FALSE)
#' mc_trends(tab, meta, "sample", "subject", "time", verbose = FALSE)
#'
#' @export
print.mc_trends <- function(x, ...) {
    mc_trend_print(x, "mc_trends")
    invisible(x)
}


#' Print a group-comparison result
#'
#' @param x An `mc_compare` object, as returned by [mc_compare()].
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#'
#' @examples
#' set.seed(1)
#' meta <- expand.grid(
#'     time = 1:6, subject = paste0("S", 1:8),
#'     stringsAsFactors = FALSE
#' )
#' meta$sample <- paste0(meta$subject, "_T", meta$time)
#' meta$arm <- ifelse(
#'     meta$subject %in% paste0("S", 1:4), "control", "treated"
#' )
#' m <- matrix(rnorm(nrow(meta)), nrow = 1,
#'             dimnames = list("Diverges", meta$sample))
#' m[1, ] <- m[1, ] + ifelse(meta$arm == "treated", 0.9 * meta$time, 0)
#' tab <- data.frame(taxon = rownames(m), m, check.names = FALSE)
#' mc_compare(
#'     tab, meta, "sample", "subject", "time",
#'     group_col = "arm", verbose = FALSE
#' )
#'
#' @export
print.mc_compare <- function(x, ...) {
    mc_trend_print(x, "mc_compare")
    invisible(x)
}


#' Shared printing for both result classes
#'
#' @param x The result object.
#' @param what The class name to show.
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
#' @noRd
mc_trend_print <- function(x, what) {
    s <- x$settings
    r <- x$results
    hit <- !is.na(r$q_value) & r$q_value < s$alpha

    cat("<", what, ">\n", sep = "")
    cat("  taxa tested  : ", sum(r$converged), " of ", nrow(r), "\n", sep = "")
    cat("  model        : s(time, k = ",
        if (length(s$k) == 1L) s$k else "chosen per taxon",
        ") + ", s$subject_effect,
        " subject effect\n",
        sep = ""
    )
    # Silent when the table was fitted as supplied, which is the default and
    # the usual case. Reported whenever it was not, since the same taxon on
    # two scales is two different results.
    if (!identical(s$normalization, "NONE") ||
        !identical(s$transform, "NONE")) {
        cat("  input scale  : ", s$normalization, " + ", s$transform, "\n",
            sep = ""
        )
    }
    if (!is.null(s$group_col)) {
        cat("  group        : ", s$group_col, " (", length(s$levels),
            " levels), reference ", s$reference, "\n",
            sep = ""
        )
        # Which test was used is decided by the number of levels, not by
        # anything the caller asked for, so it is reported rather than left
        # to be inferred.
        cat("  test         : ", s$test, "\n", sep = "")
    }
    cat("  q < ", s$alpha, "     : ", sum(hit), "\n", sep = "")
    if (any(hit)) {
        top <- utils::head(r[hit, , drop = FALSE], 5)
        cat("\n  top taxa:\n")
        for (i in seq_len(nrow(top))) {
            cat(sprintf(
                "    %-28s edf %5.2f   q %.3g\n",
                substr(top$taxon[i], 1, 28), top$edf[i], top$q_value[i]
            ))
        }
    }
    invisible(NULL)
}
