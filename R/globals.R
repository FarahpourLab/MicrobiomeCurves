#' @importFrom stats fitted quantile
NULL

# Column names used inside dplyr and ggplot2 verbs are resolved at run time
# against the data, not against the namespace. R CMD check cannot see that and
# reports each one as an undefined global variable. Declaring them here silences
# those notes. This is a declaration only: it does not affect any computation.

utils::globalVariables(c(
    ".",
    ".data",
    "imputed_value",
    "lower",
    "obs_flag",
    "outlier_flag",
    "se",
    "silhouette",
    "taxon_idx",
    "time",
    "true_value",
    "type",
    "upper",
    "value"
))
