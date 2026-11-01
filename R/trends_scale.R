############################################
# Normalisation and transformation for mc_trends() and mc_compare().
#
# The model is Gaussian with an identity link, so the scale the table
# arrives on is the caller's to declare rather than something to guess from
# the numbers. The vocabulary is MaAsLin 2's, so a result from this package
# can be set beside one from that package: normalization is one of TSS, CLR
# or NONE, transform is one of LOG, LOGIT, AST or NONE.
#
# Two axes are involved and they are not interchangeable. TSS and CLR work
# within a sample, because a proportion and a geometric mean are properties
# of one composition. LOG's zero replacement works within a taxon, because a
# detection limit belongs to the taxon. MaAsLin 2 splits them the same way:
# decostand(method = "total", MARGIN = 1) over samples in rows, then
# apply(features, 2, LOG) over features in columns. This package holds taxa
# in rows and samples in columns, so TSS and CLR run down the columns and
# LOG runs along the rows.
#
# Zero handling is deliberately not uniform, and the two rules are:
#   CLR   pseudocount, the same argument the imputation path takes, so a
#         trend test and an imputation of one table sit on the same values.
#   LOG   half the smallest positive value in the taxon, which is MaAsLin 2's
#         own rule, so TSS + LOG reproduces that package rather than
#         approximating it.
############################################


#' Reject a normalisation and transform pair that cannot both apply
#'
#' @description
#' Only one pair is impossible rather than merely unwise. CLR returns values
#' on the whole real line, and all three transforms need a bounded
#' non-negative input, so the combination has no meaning and is stopped
#' instead of producing `NaN`. MaAsLin 2 restricts LOGIT and AST the same
#' way, to normalisations that leave proportions behind.
#'
#' @param normalization,transform As in [mc_trends()], already resolved.
#'
#' @return `NULL`, invisibly. Called for the error.
#'
#' @keywords internal
#' @noRd
mc_trend_check_scale <- function(normalization, transform) {
    if (identical(normalization, "CLR") && !identical(transform, "NONE")) {
        stop(
            "normalization = \"CLR\" cannot be combined with transform = \"",
            transform, "\". Centred log ratios are negative, so LOG is ",
            "undefined and LOGIT and AST are outside their domain. ",
            "Use transform = \"NONE\".",
            call. = FALSE
        )
    }
    invisible(NULL)
}


#' Put an abundance matrix onto the scale the caller asked for
#'
#' @description
#' Applied to the whole table before any taxon is selected, because both
#' normalisations are properties of the sample as a whole: a proportion
#' needs every taxon's share and a centred log ratio needs the geometric
#' mean over the whole composition. Normalising a subset would answer a
#' different question from the one the caller asked.
#'
#' @param mat Numeric matrix, taxa in rows and samples in columns.
#' @param normalization,transform,pseudocount As in [mc_trends()].
#'
#' @return A matrix the same shape as `mat`.
#'
#' @keywords internal
#' @noRd
mc_trend_scale <- function(mat, normalization, transform, pseudocount) {
    mc_trend_check_scale(normalization, transform)

    if (identical(normalization, "TSS")) {
        mat <- mc_trend_tss(mat)
    } else if (identical(normalization, "CLR")) {
        mc_check_pseudocount(pseudocount)
        mat <- mc_clr(mat, pseudocount)
    }

    if (identical(transform, "LOG")) {
        mat <- mc_trend_log(mat)
    } else if (identical(transform, "LOGIT")) {
        mat <- mc_trend_logit(mat)
    } else if (identical(transform, "AST")) {
        mat <- mc_trend_ast(mat)
    }

    mat
}


#' How the chosen scale reads, for the console and for the log
#'
#' @description
#' One formatter with two readers, because a scale reported one way on
#' screen and another way in the file would make the two hard to compare
#' when a result is being checked months later.
#'
#' @param normalization,transform,pseudocount As in [mc_trends()].
#'
#' @return A character vector of aligned lines.
#'
#' @keywords internal
#' @noRd
mc_trend_scale_lines <- function(normalization, transform, pseudocount) {
    norm_note <- switch(normalization,
        NONE = "NONE, the table is fitted as supplied",
        TSS = "TSS, proportions within each sample",
        CLR = paste0(
            "CLR, within each sample, pseudocount ",
            if (identical(pseudocount, "auto")) {
                paste0(
                    "\"auto\" (", mc_zero_mult,
                    " of the smallest positive value)"
                )
            } else {
                paste0(pseudocount, ", added to every entry")
            }
        )
    )
    trans_note <- switch(transform,
        NONE = "NONE",
        LOG = "LOG, log2, zeros halved within each taxon",
        LOGIT = "LOGIT, with 0 and 1 set to 0 as MaAsLin 2 does",
        AST = "AST, arcsine square root"
    )

    c(
        paste0("  normalization  : ", norm_note),
        paste0("  transform      : ", trans_note)
    )
}


#' Total sum scaling, within each sample
#'
#' @param mat Numeric matrix, taxa in rows and samples in columns.
#'
#' @return A matrix of proportions, each column summing to one over the
#'   taxa that are not `NA`.
#'
#' @keywords internal
#' @noRd
mc_trend_tss <- function(mat) {
    if (any(mat < 0, na.rm = TRUE)) {
        stop(
            "normalization = \"TSS\" needs non-negative values: abundance ",
            "holds negative numbers, so it cannot be an untransformed ",
            "abundance table. If it is already on an analysis scale, use ",
            "normalization = \"NONE\".",
            call. = FALSE
        )
    }

    tot <- colSums(mat, na.rm = TRUE)
    out <- mat

    # A sample holding nothing has no total to divide by. Dividing anyway
    # gives NaN, and a GAM drops NA rows while NaN would travel further, so
    # the sample leaves as NA and the fit simply does without it.
    bad <- !is.finite(tot) | tot <= 0
    if (any(!bad)) {
        out[, !bad] <- sweep(
            mat[, !bad, drop = FALSE], 2, tot[!bad], "/"
        )
    }
    if (any(bad)) {
        out[, bad] <- NA_real_
    }
    out
}


#' Log2 with MaAsLin 2's zero replacement, within each taxon
#'
#' @description
#' Zeros are replaced by half the smallest positive value the taxon reaches
#' anywhere in the study, which is MaAsLin 2's `LOG`. The axis is the taxon
#' and not the sample, matching `apply(features, 2, LOG)` on a table of
#' samples by features, and the constant is theirs rather than the
#' package's, so that `normalization = "TSS"` with this transform can be
#' compared against a MaAsLin 2 result directly.
#'
#' @param mat Numeric matrix, taxa in rows and samples in columns.
#'
#' @return A matrix of log2 values.
#'
#' @keywords internal
#' @noRd
mc_trend_log <- function(mat) {
    if (any(mat < 0, na.rm = TRUE)) {
        stop(
            "transform = \"LOG\" needs non-negative values: abundance holds ",
            "negative numbers, so the logarithm is undefined. If the table ",
            "is already on a log or CLR scale, use transform = \"NONE\".",
            call. = FALSE
        )
    }

    out <- mat
    for (i in seq_len(nrow(mat))) {
        x <- mat[i, ]
        pos <- x[!is.na(x) & x > 0]
        if (length(pos) == 0) {
            # A taxon measured as zero everywhere has no scale of its own to
            # draw a replacement from. It also carries nothing to fit, so it
            # leaves as NA and is reported as a taxon that could not be
            # tested rather than stopping the run.
            out[i, ] <- NA_real_
            next
        }
        x[!is.na(x) & x == 0] <- min(pos) / 2
        out[i, ] <- log2(x)
    }
    out
}


#' Logit, as MaAsLin 2 applies it
#'
#' @description
#' `car::logit(x, adjust = 0)` adds nothing, so an exact zero or one gives
#' an infinite result, and MaAsLin 2 then sets every non-finite value to
#' zero. That is kept here for comparability and warned about, because a
#' logit of zero is a proportion of one half: an exact zero does not land at
#' the bottom of the scale but in the middle of it.
#'
#' @param mat Numeric matrix of proportions, taxa in rows.
#'
#' @return A matrix of logits.
#'
#' @keywords internal
#' @noRd
mc_trend_logit <- function(mat) {
    if (any(mat < 0 | mat > 1, na.rm = TRUE)) {
        stop(
            "transform = \"LOGIT\" needs proportions between 0 and 1. Use ",
            "normalization = \"TSS\" first, or transform = \"NONE\".",
            call. = FALSE
        )
    }

    out <- log(mat / (1 - mat))
    hit <- !is.na(out) & !is.finite(out)
    if (any(hit)) {
        warning(
            sum(hit), " value(s) were exactly 0 or 1, where the logit is ",
            "infinite. They are set to 0, as MaAsLin 2 does, which places ",
            "them at a proportion of one half rather than at the end of ",
            "the scale. Consider transform = \"AST\", which has no such ",
            "boundary.",
            call. = FALSE
        )
        out[hit] <- 0
    }
    out
}


#' Arcsine square root, as MaAsLin 2 applies it
#'
#' @param mat Numeric matrix, taxa in rows, values between -1 and 1.
#'
#' @return A matrix the same shape as `mat`.
#'
#' @keywords internal
#' @noRd
mc_trend_ast <- function(mat) {
    if (any(abs(mat) > 1, na.rm = TRUE)) {
        stop(
            "transform = \"AST\" is defined for values between -1 and 1. ",
            "Use normalization = \"TSS\" first, or transform = \"NONE\".",
            call. = FALSE
        )
    }
    sign(mat) * asin(sqrt(abs(mat)))
}
