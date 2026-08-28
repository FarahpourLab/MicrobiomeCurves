#' Check that an object came from mc_prepare()
#'
#' @param prep The object to check.
#'
#' @return `NULL`, invisibly. Stops when a required element is absent.
#'
#' @keywords internal
#' @noRd
mc_check_prep <- function(prep) {
    required_names <- c(
        "dat", "taxon_col", "col_map", "reps", "times", "mask_pairs"
    )
    if (!all(required_names %in% names(prep))) {
        stop("prep must be an object returned by mc_prepare()")
    }
    invisible(NULL)
}

#' Choose the number of clusters for one taxon
#'
#' @description
#' Returns `K` when the caller fixed it, and otherwise selects it by mean
#' silhouette width.
#'
#' @param scores Numeric matrix of FPCA scores, one row per subject.
#' @param taxon_name Character name of the taxon, used in the plot title.
#' @param K Integer or `NULL`. Fixed number of clusters, or `NULL` to select.
#' @param seed Integer random seed.
#'
#' @return Integer number of clusters.
#'
#' @keywords internal
#' @noRd
mc_choose_k <- function(scores, taxon_name, K, seed) {
    if (is.null(K)) {
        select_K_silhouette_plot(scores, taxon_name, seed = seed, plot = FALSE)
    } else {
        K
    }
}

#' Build the observed trajectory of every subject for one taxon
#'
#' @param df_sp Long data frame for a single taxon, with `replicate`, `time`
#'   and `value` columns.
#' @param all_reps Character vector of every subject identifier.
#'
#' @return List with `Ly` and `Lt`, each a named list holding one numeric
#'   vector per subject. Missing observations are dropped.
#'
#' @keywords internal
#' @noRd
mc_taxon_trajectories <- function(df_sp, all_reps) {
    Lt <- Ly <- vector("list", length(all_reps))
    names(Lt) <- names(Ly) <- all_reps

    for (r in all_reps) {
        tmp <- dplyr::filter(df_sp, replicate == r)
        keep <- !is.na(tmp$value)
        Lt[[r]] <- tmp$time[keep]
        Ly[[r]] <- tmp$value[keep]
    }

    list(Ly = Ly, Lt = Lt)
}

#' Flag outlying subjects for one taxon
#'
#' @param Ly,Lt Named lists of observed values and times, one entry per
#'   subject.
#' @param use_outliers Logical. When `FALSE` no subject is flagged.
#'
#' @return Named logical vector, one element per subject.
#'
#' @keywords internal
#' @noRd
mc_taxon_outliers <- function(Ly, Lt, use_outliers) {
    outliers <- if (use_outliers) {
        detect_outliers_depth(Ly, Lt, alpha = 0.05)
    } else {
        rep(FALSE, length(Ly))
    }
    names(outliers) <- names(Ly)
    outliers
}

#' Cluster the subjects of one taxon
#'
#' @description
#' Fits FPCA to the non-outlying subjects and groups their scores with
#' k-means.
#'
#' @param Ly,Lt Named lists of observed values and times, one entry per
#'   subject.
#' @param outliers Named logical vector marking outlying subjects.
#' @param K Integer or `NULL`, passed to [mc_choose_k()].
#' @param taxon_name Character name of the taxon.
#' @param seed Integer random seed.
#'
#' @return Named integer vector giving a cluster for each subject used in the
#'   fit. Every subject is placed in cluster 1 when FPCA fails or a single
#'   cluster was requested.
#'
#' @keywords internal
#' @noRd
mc_taxon_clusters <- function(Ly, Lt, outliers, K, taxon_name, seed,
                               cluster_method = "fpca") {
    Ly_use <- Ly[!outliers]
    Lt_use <- Lt[!outliers]

    if (length(Ly_use) < 2) {
        Ly_use <- Ly
        Lt_use <- Lt
    }

    single_cluster <- function() {
        cl <- rep(1, length(Ly_use))
        names(cl) <- names(Ly_use)
        cl
    }

    fp_clust <- mc_safe_fpca(Ly_use, Lt_use)
    if (is.null(fp_clust) || is.null(fp_clust$xiEst)) {
        return(single_cluster())
    }

    scores <- as.matrix(scale(fp_clust$xiEst))
    K_use <- mc_choose_k(scores, taxon_name, K, seed)

    if (K_use <= 1 || nrow(scores) < 2) {
        return(single_cluster())
    }

    mc_assign_clusters(
        fp_clust, scores, K_use, names(Ly_use), cluster_method
    )
}

#' Find the subjects sharing a cluster with the target
#'
#' @param Ly_tmp,Lt_tmp Named lists of values and times for the subjects still
#'   in play.
#' @param r_k Character identifier of the target subject.
#' @param K Integer or `NULL`, passed to [mc_choose_k()].
#' @param taxon_name Character name of the taxon.
#' @param seed Integer random seed.
#'
#' @return Character vector of subject identifiers. Every available subject is
#'   returned when clustering does not apply or the cluster would hold one
#'   member.
#'
#' @keywords internal
#' @noRd
mc_members_of_cluster <- function(Ly_tmp, Lt_tmp, r_k, K, taxon_name, seed,
                                   cluster_method = "fpca") {
    fp_tmp <- mc_safe_fpca(Ly_tmp, Lt_tmp)
    if (is.null(fp_tmp) || is.null(fp_tmp$xiEst)) {
        return(names(Ly_tmp))
    }

    scores <- as.matrix(scale(fp_tmp$xiEst))
    K_use <- mc_choose_k(scores, taxon_name, K, seed)

    if (K_use <= 1 || nrow(scores) < 2) {
        return(names(Ly_tmp))
    }

    clusters_tmp <- mc_assign_clusters(
        fp_tmp, scores, K_use, names(Ly_tmp), cluster_method
    )

    members <- if (!(r_k %in% names(clusters_tmp))) {
        names(Ly_tmp)
    } else {
        cl_id <- clusters_tmp[r_k]
        names(clusters_tmp)[clusters_tmp == cl_id]
    }

    # A cluster of one carries no information, so use every available subject.
    if (length(members) < 2) {
        members <- names(Ly_tmp)
    }

    members
}

#' Choose which subjects inform one imputed cell
#'
#' @description
#' Drops outlying subjects, keeping the target subject even when it is
#' flagged, then restricts to the target's own cluster where clustering
#' applies.
#'
#' @param Ly,Lt Named lists of observed values and times, one entry per
#'   subject.
#' @param outliers Named logical vector marking outlying subjects.
#' @param r_k Character identifier of the target subject.
#' @param K Integer or `NULL`, passed to [mc_choose_k()].
#' @param taxon_name Character name of the taxon.
#' @param seed Integer random seed.
#'
#' @return List with `members`, the subject identifiers to fit on, and
#'   `Ly_tmp` and `Lt_tmp`, the fallback set used when the cluster is too
#'   small.
#'
#' @keywords internal
#' @noRd
mc_cell_members <- function(Ly, Lt, outliers, r_k, K, taxon_name, seed,
                             cluster_method = "fpca") {
    keep_idx <- if (outliers[r_k]) {
        (!outliers) | (names(Ly) == r_k)
    } else {
        !outliers
    }

    Ly_tmp <- Ly[keep_idx]
    Lt_tmp <- Lt[keep_idx]

    if (length(Ly_tmp) < 2) {
        Ly_tmp <- Ly
        Lt_tmp <- Lt
    }

    out <- function(members) {
        list(members = members, Ly_tmp = Ly_tmp, Lt_tmp = Lt_tmp)
    }

    if (length(Ly_tmp) < 2) {
        return(out(names(Ly_tmp)))
    }

    out(mc_members_of_cluster(
        Ly_tmp, Lt_tmp, r_k, K, taxon_name, seed, cluster_method
    ))
}

#' Impute one subject-timepoint for one taxon
#'
#' @param Ly,Lt Named lists of observed values and times, one entry per
#'   subject.
#' @param sel List returned by [mc_cell_members()].
#' @param r_k Character identifier of the target subject.
#' @param tt Numeric time point to predict.
#'
#' @return List from [mc_analytic_ci()], or `NULL` when FPCA could not be
#'   fitted.
#'
#' @keywords internal
#' @noRd
mc_impute_cell <- function(Ly, Lt, sel, r_k, tt) {
    Ly_cl <- Ly[sel$members]
    Lt_cl <- Lt[sel$members]

    if (length(Ly_cl) < 2) {
        Ly_cl <- sel$Ly_tmp
        Lt_cl <- sel$Lt_tmp
    }

    fp <- mc_safe_fpca(Ly_cl, Lt_cl)
    if (is.null(fp)) {
        return(NULL)
    }

    mc_analytic_ci(
        fp,
        obs_times = Lt[[r_k]],
        obs_values = Ly[[r_k]],
        pred_time = tt,
        include_noise = TRUE
    )
}

#' Fit one taxon and impute each of its masked cells
#'
#' @param df_sp Long data frame for a single taxon.
#' @param reps Character vector of every subject identifier.
#' @param rows Integer row indices of `pred_long` belonging to this taxon.
#' @param pred_long Table of cells to impute.
#' @param use_outliers Logical. Whether to screen outlying trajectories.
#' @param K Integer or `NULL`. Number of clusters.
#' @param taxon_name Character name of the taxon.
#' @param seed Integer random seed.
#'
#' @return List with `outliers`, `clusters`, and `values`, a data frame of the
#'   imputed value and FPCA flag for each row in `rows`.
#'
#' @keywords internal
#' @noRd
mc_fit_taxon <- function(df_sp, reps, rows, pred_long, use_outliers,
                          K, taxon_name, seed, cluster_method = "fpca") {
    traj <- mc_taxon_trajectories(df_sp, reps)
    Ly <- traj$Ly
    Lt <- traj$Lt

    outliers <- mc_taxon_outliers(Ly, Lt, use_outliers)
    clusters <- mc_taxon_clusters(
        Ly, Lt, outliers, K, taxon_name, seed, cluster_method
    )

    imputed <- rep(NA_real_, length(rows))
    used <- rep(FALSE, length(rows))

    for (j in seq_along(rows)) {
        k <- rows[j]
        r_k <- pred_long$rep[k]
        tt <- pred_long$time[k]

        sel <- mc_cell_members(
            Ly, Lt, outliers, r_k, K, taxon_name, seed, cluster_method
        )
        ci_obj <- mc_impute_cell(Ly, Lt, sel, r_k, tt)

        if (!is.null(ci_obj)) {
            imputed[j] <- ci_obj$mean
            used[j] <- TRUE
        }
    }

    list(
        outliers = outliers,
        clusters = clusters,
        values = data.frame(imputed = imputed, used = used)
    )
}

#' Fit every taxon in turn
#'
#' @param long Long form of the masked table.
#' @param taxon_col Character name of the taxon column.
#' @param species_vec Character vector of taxon names, in row order.
#' @param reps Character vector of subject identifiers.
#' @param pred_long Table of cells to impute.
#' @param use_outliers Logical. Whether to screen outlying curves.
#' @param K Integer or `NULL`. Number of clusters.
#' @param seed Integer random seed.
#'
#' @return List with the filled `pred_long`, and `clusters` and `outliers`,
#'   each a list with one entry per taxon.
#'
#' @keywords internal
#' @noRd
mc_fit_all_taxa <- function(long, taxon_col, species_vec, reps,
                             pred_long, use_outliers, K, seed,
                             cluster_method = "fpca") {
    n_taxa <- length(species_vec)
    clusters_by_taxon <- vector("list", length = n_taxa)
    outliers_by_taxon <- vector("list", length = n_taxa)
    names(clusters_by_taxon) <- species_vec
    names(outliers_by_taxon) <- species_vec

    for (i in seq_len(n_taxa)) {
        sp_name <- species_vec[i]
        df_sp <- dplyr::filter(long, .data[[taxon_col]] == sp_name)
        rows <- which(pred_long$taxon_idx == i)

        res <- mc_fit_taxon(
            df_sp, reps, rows, pred_long, use_outliers, K, sp_name, seed,
            cluster_method
        )

        outliers_by_taxon[[i]] <- res$outliers
        clusters_by_taxon[[i]] <- res$clusters
        pred_long$imputed_value[rows] <- res$values$imputed
        pred_long$FPCA_used[rows] <- res$values$used
    }

    list(
        pred_long = pred_long,
        clusters = clusters_by_taxon,
        outliers = outliers_by_taxon
    )
}

#' Mask the target cells and reshape the table to long form
#'
#' @param dat Wide abundance table.
#' @param col_map Column map from [mc_prepare()].
#' @param mask_pairs Rows of `col_map` to mask.
#'
#' @return List with `dat`, the masked wide table, and `long`, its long form
#'   with `replicate`, `time` and `value` columns.
#'
#' @keywords internal
#' @noRd
mc_mask_and_reshape <- function(dat, col_map, mask_pairs) {
    for (cl in unique(mask_pairs$col)) {
        dat[[cl]] <- NA
    }

    long <- tidyr::pivot_longer(
        dat,
        cols = dplyr::all_of(col_map$col),
        names_to = "col",
        values_to = "value"
    )
    long <- dplyr::left_join(long, col_map, by = "col")
    long <- dplyr::rename(long, replicate = rep)

    list(dat = dat, long = long)
}

#' Build the table of cells to impute
#'
#' @param dat_orig Unmasked wide abundance table.
#' @param taxon_col Character name of the taxon column.
#' @param mask_pairs Masked subject-timepoint rows from [mc_prepare()].
#' @param n_taxa Integer number of taxa.
#'
#' @return Data frame with one row per taxon and masked cell, carrying the
#'   true value where one exists.
#'
#' @keywords internal
#' @noRd
mc_build_targets <- function(dat_orig, taxon_col, mask_pairs, n_taxa) {
    true_long <- tidyr::expand_grid(
        taxon_idx = seq_len(n_taxa),
        mask_pairs
    )

    dplyr::mutate(
        true_long,
        species = dat_orig[[taxon_col]][taxon_idx],
        true_value = mapply(
            function(i, col) dat_orig[[col]][i], taxon_idx, col
        ),
        imputed_value = NA_real_,
        FPCA_used = FALSE
    )
}
