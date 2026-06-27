# BAMITA dependencies. The functions file loads no packages itself.
library(rTensor)
library(tensor)
library(multiway)
library(MASS)
library(bayesm)
library(matrixNormal)
library(abind)
library(expm)
library(MCMCpack)
library(MBSP)
library(invgamma)
library(dae)
library(parallel)

set.seed(12345)

# ---------------------------------------------------------------------------
# BAMITA for the BONUS (CF infants) dataset.
# Single cohort (one sorted input, no group token), 7 timepoints. All three
# missingness mechanisms (MCAR, MAR, MNAR) are processed. No parameter tuning.
# This script lives 3 levels under the repo root.
# ---------------------------------------------------------------------------

repo_root <- normalizePath(file.path(getwd(), "..", "..", ".."), winslash = "/", mustWork = TRUE)
functions_file <- file.path(repo_root, "02.Competitors", "MultiwayImputation", "Bayesian_Tensor_Imputation_Functions.R")
source(functions_file)

dataset <- "BONUS"
input_dir <- file.path(repo_root, "01.Benchmarking Datasets", dataset, "03.processed", "04.sorted")
mechanisms <- c("MCAR", "MAR", "MNAR")
tax_level_expected <- "l7_species"

# tutorial i.i.d. settings (no tuning)
num_com <- 3
max_iter <- 5000
burn_in <- 2000
thin <- 1

parse_tax_level <- function(file_name) {
  parts <- strsplit(sub("^sorted_clr_transformation_relative_abundance_", "", file_name), "_", fixed = TRUE)[[1]]
  paste(parts[1:2], collapse = "_")
}

parse_pseudo_count <- function(file_name) {
  stem <- sub("\\.csv$", "", file_name)
  parts <- strsplit(sub("^sorted_clr_transformation_relative_abundance_", "", stem), "_", fixed = TRUE)[[1]]
  as.numeric(parts[3])   # BONUS has no group token: <tax>_<tax2>_<pseudo>
}

extract_tensor <- function(df) {
  sample_cols <- setdiff(colnames(df), "cluster")
  sample_parts <- do.call(rbind, strsplit(sample_cols, "\\.", fixed = FALSE))
  col_map <- data.frame(
    sample = sample_cols,
    subject = sample_parts[, 1],
    timepoint = as.integer(sample_parts[, 2]),
    stringsAsFactors = FALSE
  )
  col_map <- col_map[order(as.integer(sub("^[A-Za-z]+", "", col_map$subject)), col_map$timepoint), ]

  subjects <- unique(col_map$subject)
  timepoints <- sort(unique(col_map$timepoint))

  tensor <- array(
    NA_character_,
    dim = c(nrow(df), length(subjects), length(timepoints)),
    dimnames = list(rownames(df), subjects, as.character(timepoints))
  )

  for (row_idx in seq_len(nrow(col_map))) {
    subject_idx <- match(col_map$subject[row_idx], subjects)
    time_idx <- match(col_map$timepoint[row_idx], timepoints)
    tensor[, subject_idx, time_idx] <- as.character(df[[col_map$sample[row_idx]]])
  }

  list(tensor = tensor, subjects = subjects, timepoints = timepoints)
}

mask_tensor <- function(tensor, mask_df, subjects, timepoints) {
  mis <- array(0L, dim = dim(tensor), dimnames = dimnames(tensor))
  mask_subject_col <- names(mask_df)[1]
  mask_time_cols <- setdiff(names(mask_df), mask_subject_col)

  for (row_idx in seq_len(nrow(mask_df))) {
    subject <- as.character(mask_df[[mask_subject_col]][row_idx])
    subject_idx <- match(subject, subjects)
    if (is.na(subject_idx)) next

    for (time_name in mask_time_cols) {
      time_idx <- match(as.integer(time_name), timepoints)
      if (is.na(time_idx)) next

      if (as.numeric(mask_df[[time_name]][row_idx]) == 0) {
        tensor[, subject_idx, time_idx] <- "NA"
        mis[, subject_idx, time_idx] <- 1L
      }
    }
  }

  list(tensor = tensor, mis = mis)
}

posterior_mean_tensor <- function(fit, mis) {
  sample_count <- dim(fit$resu1)[1]
  tensor_shape <- c(dim(fit$resu1)[2], dim(fit$resu2)[2], dim(fit$resu3)[2])
  posterior_samples <- array(NA_real_, dim = c(sample_count, tensor_shape))

  for (idx in seq_len(sample_count)) {
    u_est <- list(fit$resu1[idx, , ], fit$resu2[idx, , ], fit$resu3[idx, , ])
    sample_tensor <- outerimp(u_est)
    sample_tensor <- sample_tensor + rnorm(prod(tensor_shape), mean = 0, sd = sqrt(fit$sigma_err[idx])) * mis
    posterior_samples[idx, , , ] <- sample_tensor
  }

  apply(posterior_samples, c(2, 3, 4), mean)
}

tensor_to_wide <- function(tensor, sample_order) {
  out <- matrix(NA_real_, nrow = dim(tensor)[1], ncol = length(sample_order))
  colnames(out) <- sample_order

  split_names <- do.call(rbind, strsplit(sample_order, "\\.", fixed = FALSE))
  subjects <- unique(split_names[, 1])
  timepoints <- sort(unique(as.integer(split_names[, 2])))

  for (sample_name in sample_order) {
    bits <- strsplit(sample_name, "\\.", fixed = FALSE)[[1]]
    subject_idx <- match(bits[1], subjects)
    time_idx <- match(as.integer(bits[2]), timepoints)
    out[, sample_name] <- tensor[, subject_idx, time_idx]
  }

  out
}

metric_summary <- function(truth_tensor, imputed_tensor, mis) {
  masked_truth <- truth_tensor[mis == 1L]
  masked_imputed <- imputed_tensor[mis == 1L]
  diffs <- masked_truth - masked_imputed

  data.frame(
    mse = mean(diffs ^ 2),
    rmse = sqrt(mean(diffs ^ 2)),
    mae = mean(abs(diffs)),
    n_missing = length(diffs),
    stringsAsFactors = FALSE
  )
}

max_cores <- 20L
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_cores) || n_cores < 1L) {
  n_cores <- max(1L, parallel::detectCores())
}
n_cores <- min(n_cores, max_cores)
cat(sprintf("Using %d core(s) for parallel mask-file processing.\n", n_cores))

RNGkind("L'Ecuyer-CMRG")
set.seed(12345)

process_mask <- function(mask_file, truth_tensor, tensor_info, truth_df,
                         file_name, tax_level, pseudo_count, mechanism, output_dir) {
  tryCatch({
    mask_ratio <- basename(dirname(mask_file))
    mask_df <- read.csv(mask_file, check.names = FALSE)

    masked <- mask_tensor(truth_tensor, mask_df, tensor_info$subjects, tensor_info$timepoints)
    observed_tensor <- masked$tensor
    mis <- masked$mis

    fit <- cpbayeimp_Jef_check_effi(
      observed_tensor, num_com = num_com, max = max_iter, burn.in = burn_in, thin = thin
    )

    imputed_tensor <- posterior_mean_tensor(fit, mis)
    output_tensor <- imputed_tensor
    observed_cells <- mis == 0L
    output_tensor[observed_cells] <- as.numeric(truth_tensor)[observed_cells]

    imputed_wide <- tensor_to_wide(output_tensor, colnames(truth_df))
    imputed_df <- data.frame(
      `#OTU ID` = rownames(truth_df),
      as.data.frame(imputed_wide, check.names = FALSE),
      check.names = FALSE
    )

    sel_out_dir <- file.path(output_dir, tax_level, mask_ratio)
    dir.create(sel_out_dir, recursive = TRUE, showWarnings = FALSE)

    out_stub <- sub("\\.csv$", "", file_name)
    mask_base <- sub("\\.csv$", "", basename(mask_file))
    mask_short <- paste(strsplit(mask_base, ".", fixed = TRUE)[[1]][-1], collapse = ".")
    out_file <- file.path(sel_out_dir, paste0("imputed_", out_stub, "_", mask_short, ".csv"))
    metrics_file <- file.path(sel_out_dir, paste0("metrics_", out_stub, "_", mask_short, ".csv"))

    write.csv(imputed_df, out_file, row.names = FALSE)

    metrics_df <- metric_summary(as.numeric(truth_tensor), as.numeric(imputed_tensor), mis)
    metrics_df$mechanism <- mechanism
    metrics_df$tax_level <- tax_level
    metrics_df$mask_ratio <- mask_ratio
    metrics_df$pseudo_count <- pseudo_count
    metrics_df$num_com <- num_com
    metrics_df$max_iter <- max_iter
    metrics_df$burn_in <- burn_in
    metrics_df$source_file <- file_name
    metrics_df$mask_file <- basename(mask_file)
    write.csv(metrics_df, metrics_file, row.names = FALSE)

    metrics_file
  }, error = function(e) {
    structure(conditionMessage(e), class = "try-error", mask_file = mask_file)
  })
}

# the single sorted input file for this dataset
input_file <- list.files(
  input_dir,
  pattern = paste0("^sorted_clr_transformation_relative_abundance_", tax_level_expected, "_.*\\.csv$"),
  full.names = TRUE
)
if (length(input_file) == 0) stop(sprintf("no sorted input found in %s", input_dir))
input_file <- input_file[1]

file_name <- basename(input_file)
tax_level <- parse_tax_level(file_name)
pseudo_count <- parse_pseudo_count(file_name)

truth_df <- read.csv(input_file, row.names = 1, check.names = FALSE)
truth_df$cluster <- NULL
tensor_info <- extract_tensor(truth_df)
truth_tensor <- tensor_info$tensor

for (mechanism in mechanisms) {
  mask_root <- file.path(repo_root, "01.Benchmarking Datasets", dataset, "03.processed", mechanism, "mask")
  output_dir <- file.path(repo_root, "04.Benchmarking", dataset, "MultiwayImputation", mechanism)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  mask_dir <- file.path(mask_root, tax_level)
  mask_files <- list.files(
    mask_dir,
    pattern = paste0("masked_metadata_", tax_level, "\\.rep_.*\\.csv$"),
    full.names = TRUE, recursive = TRUE
  )

  cat(sprintf("[%s] %s: %d mask files across %d core(s)\n",
              mechanism, file_name, length(mask_files), n_cores))

  results <- parallel::mclapply(
    mask_files, process_mask,
    truth_tensor = truth_tensor, tensor_info = tensor_info, truth_df = truth_df,
    file_name = file_name, tax_level = tax_level, pseudo_count = pseudo_count,
    mechanism = mechanism, output_dir = output_dir,
    mc.cores = n_cores, mc.preschedule = FALSE
  )

  failed <- Filter(function(x) inherits(x, "try-error"), results)
  if (length(failed) > 0) {
    warning(sprintf("[%s] %d of %d mask fits failed", mechanism, length(failed), length(mask_files)))
    for (f in failed) cat("  FAILED:", attr(f, "mask_file"), "->", as.character(f), "\n")
  }
}
