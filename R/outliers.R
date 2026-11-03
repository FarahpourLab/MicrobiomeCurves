# Outlier screening as a step of its own.
#
# The imputation screens trajectories before it clusters them, and leaves
# the flags in fit$outliers. That screening exists to serve the
# imputation: it runs on the masked table, it reports nothing, and it
# cannot be asked for without imputing something.
#
# mc_detect_outliers() asks the same question as an analysis. It returns
# every subject's depth and not only the flag, because the flag is a rank
# and not a test: the cutoff is the alpha quantile of the depths, so some
# fraction of subjects is always flagged, whether or not any trajectory is
# genuinely aberrant. A subject at the edge of a cohort is usually a real
# phenotype rather than a bad measurement, and the depth is what lets a
# reader judge which they are looking at.


#' Depth of every trajectory
#'
#' @description
#' Fraiman-Muniz depth of each curve among the others, from the FPCA fit
#' of the set. Depth measures how central a curve is: a curve near the
#' middle of the bundle scores high, one at its edge scores low.
#'
#' @param Ly List of numeric vectors, one per subject.
#' @param Lt List of numeric vectors of matching observation times.
#'
#' @return Numeric vector, one element per element of `Ly`, named as `Ly`
#'   is. A subject FPCA could not use is `NA`: it has no depth rather than
#'   a low one.
#'
#' @keywords internal
#' @noRd
mc_depth_scores <- function(Ly, Lt) {
    if (!requireNamespace("fda.usc", quietly = TRUE)) {
        stop("Package 'fda.usc' required for depth-based outliers")
    }

    out <- rep(NA_real_, length(Ly))
    names(out) <- names(Ly)

    fp <- mc_safe_fpca(Ly, Lt)
    if (is.null(fp)) {
        return(out)
    }

    curves <- stats::fitted(fp)
    grid <- fp$workGrid
    if (is.null(curves) || !is.matrix(curves) || is.null(grid) ||
        ncol(curves) != length(grid)) {
        return(out)
    }

    depth_vals <- fda.usc::depth.FM(
        fda.usc::fdata(curves, argvals = grid)
    )$dep

    # mc_safe_fpca() fits only the subjects with at least two observations,
    # so the depths are indexed against those and not against every subject
    # in Ly. Map them back to full length. When every subject was usable
    # this is the identity, which is the case the benchmark drivers always
    # produce.
    if (length(depth_vals) == length(Ly)) {
        out[] <- depth_vals
        return(out)
    }

    obs_count <- vapply(Ly, function(x) sum(!is.na(x)), numeric(1))
    keep <- which(obs_count >= 2)
    if (length(depth_vals) == length(keep)) {
        out[keep] <- depth_vals
    }
    out
}


#' Find trajectories that sit apart from the rest
#'
#' @description
#' Scores every subject's trajectory by functional depth, within each
#' taxon, and flags the least central ones. Nothing is imputed and nothing
#' is removed: the result is a report.
#'
#' @details
#' Depth is Fraiman-Muniz depth of the FPCA-fitted curve among the others
#' of that taxon. A curve near the middle of the bundle scores high and
#' one at the edge scores low, so the ordering is the finding and the flag
#' is a convenience.
#'
#' The flag marks the subjects below the `alpha` quantile of the depths.
#' That is a rank and not a test: roughly `alpha` of the subjects are
#' flagged in every taxon, whether or not any trajectory is genuinely
#' aberrant, and a taxon where every subject behaves alike still yields
#' its quota. Read `depth` before acting on `outlier`, and treat a flagged
#' subject as one worth looking at rather than one to drop. This is the
#' same screening [mc_run()] applies when `use_outliers = TRUE`, at the
#' same default `alpha`.
#'
#' @param abundance Wide abundance table: taxa in rows, samples in
#'   columns, as [mc_trends()] takes it.
#' @param metadata Data frame describing the samples.
#' @param sample_col,subject_col,time_col Column names in `metadata`.
#' @param taxon_col The column of `abundance` holding taxon names, or
#'   `NULL` to use row names.
#' @param taxa Character vector of taxa to screen, or `NULL` for all.
#' @param subjects Character vector of subjects to use, or `NULL` for all.
#'   Depth is relative to the subjects in play, so narrowing this changes
#'   the depths of those that remain.
#' @param alpha Depth quantile below which a trajectory is flagged.
#' @param normalization,transform,pseudocount The input scale, as in
#'   [mc_trends()].
#' @param seed Random seed. The caller's own random stream is restored on
#'   return.
#' @param out_dir Directory to write `outliers.tsv` into, or `NULL`.
#' @param verbose Logical. Whether to report what was done.
#'
#' @return An object of class `mc_outliers` with
#'   \itemize{
#'     \item `results`: one row per subject and taxon, with that
#'       trajectory's `depth` and whether it was flagged.
#'     \item `n_outliers`: one row per taxon, with how many subjects it
#'       held and how many were flagged.
#'     \item `settings`: what was asked for and what was done.
#'     \item `files`: anything written.
#'   }
#'
#' @examples
#' demo <- mc_demo_data(n_taxa = 3, n_subjects = 10)
#' od <- mc_detect_outliers(
#'     demo$abundance, demo$metadata,
#'     sample_col = "sample", subject_col = "subject", time_col = "time",
#'     verbose = FALSE
#' )
#' od
#'
#' # the least central subjects first, which is the finding
#' head(od$results[order(od$results$depth), ])
#'
#' @export
mc_detect_outliers <- function(abundance, metadata,
                               sample_col, subject_col, time_col,
                               taxon_col = NULL,
                               taxa = NULL,
                               subjects = NULL,
                               alpha = 0.05,
                               normalization = c("NONE", "TSS", "CLR"),
                               transform = c("NONE", "LOG", "LOGIT", "AST"),
                               pseudocount = "auto",
                               seed = 123,
                               out_dir = NULL,
                               verbose = TRUE) {
    normalization <- match.arg(normalization)
    transform <- match.arg(transform)
    say <- function(...) if (isTRUE(verbose)) message(...)

    if (length(alpha) != 1L || !is.numeric(alpha) || is.na(alpha) ||
        alpha <= 0 || alpha >= 1) {
        stop(
            "alpha must be a single number greater than 0 and less than 1.",
            call. = FALSE
        )
    }

    design <- mc_trend_design(metadata, sample_col, subject_col, time_col)
    tm <- mc_trend_matrix(
        abundance, taxon_col, taxa, normalization, transform, pseudocount
    )
    al <- mc_trend_align(design, tm$mat, verbose)
    design <- al$design
    mat <- al$mat

    if (!is.null(subjects)) {
        gone <- setdiff(subjects, design$subject)
        if (length(gone) > 0) {
            stop(
                "No sample for subject(s) ", mc_fmt_some(gone), ".",
                call. = FALSE
            )
        }
        keep <- design$subject %in% subjects
        design <- design[keep, , drop = FALSE]
        mat <- mat[, keep, drop = FALSE]
    }

    mc_preserve_rng()
    mc_reset_rng(seed)

    say(
        "mc_detect_outliers: ", nrow(mat), " taxa, ", ncol(mat),
        " samples, ", length(unique(design$subject)), " subjects, ",
        length(unique(design$time)), " time points"
    )
    say("  depth          : Fraiman-Muniz, within each taxon")
    say("  alpha          : ", alpha, " (the least central ",
        round(100 * alpha), "% of each taxon)")
    if (!identical(normalization, "NONE") || !identical(transform, "NONE")) {
        for (ln in mc_trend_scale_lines(
            normalization, transform, pseudocount
        )) {
            say(ln)
        }
    }

    # The same curve sets the clustering command builds, since the
    # question is asked of the same curves.
    units <- mc_cluster_units(design, mat, "subjects")
    res <- do.call(rbind, lapply(names(units), function(tx) {
        cv <- units[[tx]]
        if (length(cv$Ly) == 0L) {
            return(NULL)
        }
        d <- mc_depth_scores(cv$Ly, cv$Lt)
        flag <- detect_outliers_depth(cv$Ly, cv$Lt, alpha = alpha)
        data.frame(
            taxon = tx, subject = names(cv$Ly),
            depth = as.numeric(d), outlier = flag,
            stringsAsFactors = FALSE
        )
    }))

    if (is.null(res) || nrow(res) == 0L) {
        stop(
            "Nothing to screen: no taxon held a usable trajectory.",
            call. = FALSE
        )
    }
    rownames(res) <- NULL

    sp <- split(res$outlier, res$taxon)
    n_out <- data.frame(
        taxon = names(sp),
        n_subjects = vapply(sp, length, integer(1)),
        n_outliers = vapply(sp, function(x) sum(x, na.rm = TRUE), integer(1)),
        stringsAsFactors = FALSE
    )
    rownames(n_out) <- NULL

    say(
        "  screened       : ", nrow(res), " trajectory(ies) in ",
        nrow(n_out), " taxa"
    )
    say("  flagged        : ", sum(res$outlier, na.rm = TRUE))

    files <- character(0)
    if (!is.null(out_dir)) {
        files <- mc_trend_write(res, out_dir, "outliers")
        say("  written        : ", files)
    }

    structure(
        list(
            results = res,
            n_outliers = n_out,
            design = design,
            abundance = mat,
            settings = list(
                alpha = alpha, seed = seed,
                normalization = normalization, transform = transform,
                pseudocount = pseudocount,
                sample_col = sample_col, subject_col = subject_col,
                time_col = time_col,
                depth = "Fraiman-Muniz, on the FPCA fit of each taxon",
                # Said here because the flag reads like a verdict and is
                # not one: the cutoff is a quantile of the depths.
                note = paste0(
                    "the flag is the lowest ", round(100 * alpha),
                    "% by depth within each taxon, not a test"
                )
            ),
            files = files
        ),
        class = "mc_outliers"
    )
}


#' Print an outlier screening
#'
#' @param x An `mc_outliers` object, as returned by [mc_detect_outliers()].
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#'
#' @examples
#' demo <- mc_demo_data(n_taxa = 3, n_subjects = 10)
#' od <- mc_detect_outliers(
#'     demo$abundance, demo$metadata,
#'     sample_col = "sample", subject_col = "subject", time_col = "time",
#'     verbose = FALSE
#' )
#' od
#'
#' @export
print.mc_outliers <- function(x, ...) {
    s <- x$settings
    r <- x$results

    cat("<mc_outliers>\n")
    cat("  depth        : ", s$depth, "\n", sep = "")
    cat("  alpha        : ", s$alpha, "\n", sep = "")
    if (!identical(s$normalization, "NONE") ||
        !identical(s$transform, "NONE")) {
        cat("  input scale  : ", s$normalization, " + ", s$transform, "\n",
            sep = ""
        )
    }
    cat("  screened     : ", nrow(r), " trajectory(ies) in ",
        nrow(x$n_outliers), " taxa\n",
        sep = ""
    )
    cat("  flagged      : ", sum(r$outlier, na.rm = TRUE), "\n", sep = "")
    cat("  note         : ", s$note, "\n", sep = "")
    invisible(x)
}
