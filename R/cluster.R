# Clustering as a step of its own.
#
# The imputation clusters subjects inside every taxon before it fills a
# cell, and leaves the result in fit$clusters. That clustering exists to
# serve the imputation: it is computed on the masked table, it is not
# reported, and it cannot be asked for without imputing something.
#
# mc_cluster() asks the same question as an analysis rather than as a
# step. It takes the table as it stands, groups whichever curves were
# asked for, and reports the assignment. Three things can be grouped and
# they are different questions, not three spellings of one:
#
#   by = "subjects"   within one taxon, which subjects move alike
#   by = "taxa"       within one subject, which taxa move alike
#   by = "taxa_mean"  across the cohort, which taxa move alike on average
#
# All three reduce to the same operation on a set of curves, so all three
# go through the machinery the imputation already uses, and a clustering
# asked for here is the clustering that would have been used there.


#' Group longitudinal trajectories
#'
#' @description
#' Clusters trajectories without imputing anything. The curves to group
#' are chosen by `by`: subjects within a taxon, taxa within a subject, or
#' taxa by their cohort-average trajectory.
#'
#' @details
#' Each set of curves is fitted by FPCA and grouped by k-means, which is
#' the route the imputation uses internally. `method = "fpca"` clusters
#' the standardised FPC scores, so curves are compared through the leading
#' modes of variation. `method = "kmeans_fd"` clusters the fitted curves
#' themselves under an L2 metric; it needs `fda.usc`, and it falls back to
#' the score-based route on sets of curves too small for it rather than
#' losing the unit.
#'
#' `C = NULL` chooses the number of clusters by mean silhouette width,
#' separately for every unit, because the number of groups a taxon falls
#' into is a property of that taxon and not of the run. A curve set that
#' FPCA cannot fit, or that holds fewer than two curves, is reported as a
#' single cluster rather than dropped, so every unit asked for appears in
#' the result.
#'
#' @param abundance Wide abundance table: taxa in rows, samples in
#'   columns, as [mc_trends()] takes it.
#' @param metadata Data frame describing the samples.
#' @param sample_col,subject_col,time_col Column names in `metadata`.
#' @param taxon_col The column of `abundance` holding taxon names, or
#'   `NULL` to use row names.
#' @param by What to group. `"subjects"` clusters subjects within each
#'   taxon, `"taxa"` clusters taxa within each subject, and `"taxa_mean"`
#'   clusters taxa by their mean trajectory over subjects.
#' @param taxa Character vector of taxa to use, or `NULL` for all. With
#'   `by = "subjects"` this is also the set of clusterings produced.
#' @param subjects Character vector of subjects to use, or `NULL` for all.
#'   With `by = "taxa"` this is also the set of clusterings produced.
#' @param C Number of clusters, or `NULL` to choose one per unit by
#'   silhouette width.
#' @param method `"fpca"` or `"kmeans_fd"`.
#' @param normalization,transform,pseudocount The input scale, as in
#'   [mc_trends()].
#' @param seed Random seed. The caller's own random stream is restored on
#'   return.
#' @param out_dir Directory to write `clusters.tsv` into, or `NULL`.
#' @param verbose Logical. Whether to report what was done.
#'
#' @return An object of class `mc_cluster` with
#'   \itemize{
#'     \item `results`: one row per clustered curve, holding the unit it
#'       was clustered in, the curve's name and its cluster.
#'     \item `n_clusters`: one row per unit, with how many clusters it
#'       fell into and how many curves it held.
#'     \item `settings`: what was asked for and what was done.
#'     \item `files`: anything written.
#'   }
#'
#' @examples
#' demo <- mc_demo_data(n_taxa = 3, n_subjects = 8)
#'
#' # which subjects move alike, within each taxon
#' cl <- mc_cluster(
#'     demo$abundance, demo$metadata,
#'     sample_col = "sample", subject_col = "subject", time_col = "time",
#'     verbose = FALSE
#' )
#' head(cl$results)
#'
#' # which taxa move alike, on cohort average
#' ct <- mc_cluster(
#'     demo$abundance, demo$metadata,
#'     sample_col = "sample", subject_col = "subject", time_col = "time",
#'     by = "taxa_mean", verbose = FALSE
#' )
#' ct$results
#'
#' @export
mc_cluster <- function(abundance, metadata,
                       sample_col, subject_col, time_col,
                       taxon_col = NULL,
                       by = c("subjects", "taxa", "taxa_mean"),
                       taxa = NULL,
                       subjects = NULL,
                       C = NULL,
                       method = c("fpca", "kmeans_fd"),
                       normalization = c("NONE", "TSS", "CLR"),
                       transform = c("NONE", "LOG", "LOGIT", "AST"),
                       pseudocount = "auto",
                       seed = 123,
                       out_dir = NULL,
                       verbose = TRUE) {
    by <- match.arg(by)
    method <- match.arg(method)
    normalization <- match.arg(normalization)
    transform <- match.arg(transform)
    say <- function(...) if (isTRUE(verbose)) message(...)

    if (!is.null(C) && (length(C) != 1L || !is.finite(C) || C < 1L)) {
        stop("C must be NULL or a single number of at least 1.", call. = FALSE)
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

    # Reproducible for a given seed without standing on the caller's own
    # random stream, as mc_fit() does. See R/rng.R.
    mc_preserve_rng()
    mc_reset_rng(seed)

    say(
        "mc_cluster: ", nrow(mat), " taxa, ", ncol(mat), " samples, ",
        length(unique(design$subject)), " subjects, ",
        length(unique(design$time)), " time points"
    )
    say("  grouping       : ", switch(by,
        subjects = "subjects, within each taxon",
        taxa = "taxa, within each subject",
        taxa_mean = "taxa, by mean trajectory over subjects"
    ))
    say("  method         : ", method)
    say("  C              : ", if (is.null(C)) {
        "chosen per unit by silhouette width"
    } else {
        C
    })
    if (!identical(normalization, "NONE") || !identical(transform, "NONE")) {
        for (ln in mc_trend_scale_lines(
            normalization, transform, pseudocount
        )) {
            say(ln)
        }
    }

    units <- mc_cluster_units(design, mat, by)
    res <- do.call(rbind, lapply(names(units), function(u) {
        cv <- units[[u]]
        if (length(cv$Ly) == 0L) {
            return(NULL)
        }
        cl <- mc_taxon_clusters(
            cv$Ly, cv$Lt,
            outliers = rep(FALSE, length(cv$Ly)),
            C = C, taxon_name = u, seed = seed, cluster_method = method
        )
        data.frame(
            unit = u, curve = names(cl), cluster = as.integer(cl),
            stringsAsFactors = FALSE
        )
    }))

    if (is.null(res) || nrow(res) == 0L) {
        stop(
            "Nothing to cluster: no unit held two usable curves.",
            call. = FALSE
        )
    }
    rownames(res) <- NULL

    # The generic unit/curve pair is what the machinery works in, but it
    # is not what the caller asked about, so the columns are named for the
    # question rather than for the implementation.
    names(res)[1:2] <- switch(by,
        subjects = c("taxon", "subject"),
        taxa = c("subject", "taxon"),
        taxa_mean = c("set", "taxon")
    )
    if (identical(by, "taxa_mean")) {
        res <- res[, c("taxon", "cluster"), drop = FALSE]
    }

    n_cl <- mc_cluster_counts(res, by)
    say(
        "  clustered      : ", nrow(res), " curve(s) in ", nrow(n_cl),
        " unit(s)"
    )
    say(
        "  clusters       : ", min(n_cl$n_clusters), " to ",
        max(n_cl$n_clusters), " per unit"
    )

    files <- character(0)
    if (!is.null(out_dir)) {
        files <- mc_trend_write(res, out_dir, "clusters")
        say("  written        : ", files)
    }

    structure(
        list(
            results = res,
            n_clusters = n_cl,
            design = design,
            abundance = mat,
            settings = list(
                by = by, method = method, C = C, seed = seed,
                normalization = normalization, transform = transform,
                pseudocount = pseudocount,
                sample_col = sample_col, subject_col = subject_col,
                time_col = time_col,
                # What a cluster id means depends on what was clustered,
                # and ids are not comparable between units: cluster 1 of
                # one taxon has nothing to do with cluster 1 of another.
                note = "cluster ids are local to each unit"
            ),
            files = files
        ),
        class = "mc_cluster"
    )
}


#' Build the curve sets that each clustering runs on
#'
#' @description
#' One entry per unit, each holding the named lists of values and times
#' that `mc_taxon_clusters()` takes. A curve with fewer than two points
#' carries no shape and is dropped before the fit, which is why a unit can
#' come back holding fewer curves than it has members.
#'
#' @param design The aligned sample table.
#' @param mat The aligned abundance matrix.
#' @param by One of `"subjects"`, `"taxa"`, `"taxa_mean"`.
#'
#' @return A named list of `list(Ly, Lt)`.
#'
#' @keywords internal
#' @noRd
mc_cluster_units <- function(design, mat, by) {
    keep_usable <- function(cur) {
        ok <- vapply(cur, function(x) length(x$v) >= 2L, logical(1))
        cur <- cur[ok]
        list(
            Ly = lapply(cur, function(x) x$v),
            Lt = lapply(cur, function(x) x$t)
        )
    }
    one_curve <- function(v, tt) {
        ok <- !is.na(v)
        o <- order(tt[ok])
        list(v = as.numeric(v[ok][o]), t = as.numeric(tt[ok][o]))
    }

    if (identical(by, "subjects")) {
        subs <- unique(design$subject)
        out <- lapply(rownames(mat), function(tx) {
            cur <- lapply(subs, function(s) {
                i <- which(design$subject == s)
                one_curve(mat[tx, i], design$time[i])
            })
            names(cur) <- subs
            keep_usable(cur)
        })
        names(out) <- rownames(mat)
        return(out)
    }

    if (identical(by, "taxa")) {
        subs <- unique(design$subject)
        out <- lapply(subs, function(s) {
            i <- which(design$subject == s)
            cur <- lapply(rownames(mat), function(tx) {
                one_curve(mat[tx, i], design$time[i])
            })
            names(cur) <- rownames(mat)
            keep_usable(cur)
        })
        names(out) <- subs
        return(out)
    }

    # taxa_mean: one curve per taxon, averaged over whichever subjects
    # were seen at each time. The average is over the subjects present and
    # not over a fixed denominator, so a time point nobody was sampled at
    # is absent rather than zero.
    tt <- sort(unique(design$time))
    cur <- lapply(rownames(mat), function(tx) {
        v <- vapply(tt, function(u) {
            x <- mat[tx, design$time == u]
            if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
        }, numeric(1))
        one_curve(v, tt)
    })
    names(cur) <- rownames(mat)
    list(all = keep_usable(cur))
}


#' How many clusters each unit fell into
#'
#' @param res The assignment table.
#' @param by One of `"subjects"`, `"taxa"`, `"taxa_mean"`.
#'
#' @return A data.frame of one row per unit.
#'
#' @keywords internal
#' @noRd
mc_cluster_counts <- function(res, by) {
    unit_col <- switch(by,
        subjects = "taxon", taxa = "subject", taxa_mean = NA_character_
    )
    if (is.na(unit_col)) {
        return(data.frame(
            n_curves = nrow(res),
            n_clusters = length(unique(res$cluster)),
            stringsAsFactors = FALSE
        ))
    }
    sp <- split(res$cluster, res[[unit_col]])
    out <- data.frame(
        unit = names(sp),
        n_curves = vapply(sp, length, integer(1)),
        n_clusters = vapply(sp, function(x) length(unique(x)), integer(1)),
        stringsAsFactors = FALSE
    )
    names(out)[1] <- unit_col
    rownames(out) <- NULL
    out
}


#' Print a clustering
#'
#' @param x An `mc_cluster` object, as returned by [mc_cluster()].
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#'
#' @examples
#' demo <- mc_demo_data(n_taxa = 3, n_subjects = 8)
#' cl <- mc_cluster(
#'     demo$abundance, demo$metadata,
#'     sample_col = "sample", subject_col = "subject", time_col = "time",
#'     verbose = FALSE
#' )
#' cl
#'
#' @export
print.mc_cluster <- function(x, ...) {
    s <- x$settings
    n <- x$n_clusters

    cat("<mc_cluster>\n")
    cat("  grouping     : ", switch(s$by,
        subjects = "subjects, within each taxon",
        taxa = "taxa, within each subject",
        taxa_mean = "taxa, by mean trajectory"
    ), "\n", sep = "")
    cat("  method       : ", s$method, "\n", sep = "")
    cat("  C            : ", if (is.null(s$C)) "silhouette" else s$C,
        "\n",
        sep = ""
    )
    if (!identical(s$normalization, "NONE") ||
        !identical(s$transform, "NONE")) {
        cat("  input scale  : ", s$normalization, " + ", s$transform, "\n",
            sep = ""
        )
    }
    cat("  curves       : ", nrow(x$results), " in ", nrow(n), " unit(s)\n",
        sep = ""
    )
    cat("  clusters     : ", min(n$n_clusters), " to ", max(n$n_clusters),
        " per unit\n",
        sep = ""
    )
    cat("  note         : ", s$note, "\n", sep = "")
    invisible(x)
}
