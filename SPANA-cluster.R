# SPANA-cluster: Saturating Perturbation Analysis of Noncoding Activity
#
# Copyright (C) 2026 The Regents of the University of Michigan
# Emmer Lab
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 2 of the License, or (at your option)
# any later version. It is distributed WITHOUT ANY WARRANTY; see the LICENSE
# file for details.

VERSION <- "1.0.0"  # script version; recorded in the run parameter log

###############################################################################
# SPANA-cluster — nearest-neighbor clustering test of guide hits
#
# DESCRIPTION:
#   Tests whether significant guides (hits) lie closer together along the
#   screened locus than expected by chance. The test statistic is the median
#   distance from each hit to its nearest neighboring hit. The null
#   distribution is built by repeatedly drawing the same number of guides at
#   random from all tested guides, keeping their actual positions, and
#   recomputing the statistic. The one-sided empirical p-value is
#   (number of permutations with a median distance <= observed + 1) /
#   (number of permutations + 1). Depleted and enriched hits are tested
#   separately.
#
# REQUIREMENTS:
#   Base R only.
#
# EXPECTED INPUT (DESeq2_results_test.tsv from SPANA-null):
#   Required columns:
#     GuideID        — unique guide identifier
#     Coordinate     — genomic coordinate of the guide cut site
#     type           — guide class; only rows where type == "test" are used
#     adj-p deplete  — BH-adjusted p-value for depletion
#     adj-p enrich   — BH-adjusted p-value for enrichment
#   A guide is a hit if its adjusted p-value is below the FDR threshold.
#   Guides with a missing coordinate or adjusted p-value are excluded.
#
# REPRODUCIBILITY:
#   Guides are sorted by coordinate (guide ID breaks ties) before permuting,
#   and each direction starts from RANDOM_SEED, so results are identical for
#   a given input and seed regardless of input row order.
#
# OUTPUT (written to a timestamped subdirectory):
#   - <prefix>_nn_null_distributions.tsv : permutation null values of the
#       median nearest-neighbor distance (bp), one column per direction
#   - <prefix>_nn_summary.tsv            : per direction, the number of hits
#       and guides tested, observed median nearest-neighbor distance (bp),
#       empirical p-value, and number of permutations
#   - <prefix>_run_parameters.txt        : parameters, guide counts, results,
#       and R session information
#
# WORKING DIRECTORY:
#   Set your working directory to the SPANA-null output folder before
#   sourcing this script, using either the RStudio Files panel (More -> Set As
#   Working Directory) or setwd() in the console. The script will use getwd()
#   if WORK_DIR is left as "".
#
#   To bypass all readline() prompts, pre-fill the PARAMETERS block below
#   and set INTERACTIVE_MODE <- FALSE.
###############################################################################


# ============================================================================
# PARAMETERS
# Pre-fill these and set INTERACTIVE_MODE <- FALSE to run non-interactively
# (e.g. on an HPC cluster or within a pipeline).
# ============================================================================

INTERACTIVE_MODE <- TRUE   # Set FALSE to skip all prompts

# --- Paths ---
WORK_DIR    <- ""          # Working directory; leave "" to use getwd()
RES_FILE    <- ""          # Input TSV filename (leave "" to auto-detect)
OUT_PREFIX  <- ""          # Output file prefix (e.g. "MyScreen")

# --- Analysis parameters ---
N_PERM      <- 100000L     # Number of permutations per direction
FDR_THR     <- 0.01        # Hits: BH-adjusted p-value below this threshold
RANDOM_SEED <- 1L          # Random seed for reproducibility

# ============================================================================
# SETUP
# ============================================================================

if (nchar(trimws(WORK_DIR)) > 0) {
  setwd(WORK_DIR)
}
WORK_DIR <- getwd()
cat("Working directory:\n ", getwd(), "\n")

if (INTERACTIVE_MODE) {
  cat("\nDirectory contents:\n")
  print(list.files(WORK_DIR))
  if (toupper(readline("Continue? (Y/N): ")) != "Y") stop("Stopped by user.")
}

# ============================================================================
# OUTPUT PREFIX AND INTERACTIVE OVERRIDES
# ============================================================================

if (INTERACTIVE_MODE && nchar(trimws(OUT_PREFIX)) == 0) {
  OUT_PREFIX <- readline("Output prefix: ")
}
if (nchar(trimws(OUT_PREFIX)) == 0) stop("Output prefix cannot be empty.")

now     <- Sys.time()
ts      <- format(now, "%Y_%m%d_%H%M")
out_dir <- paste0(ts, "_cluster_output")

if (INTERACTIVE_MODE) {

  # FDR threshold
  repeat {
    fdr_in <- readline(paste0("FDR threshold for hits (default ", FDR_THR, "): "))
    if (nchar(trimws(fdr_in)) == 0) break
    tmp <- suppressWarnings(as.numeric(fdr_in))
    if (!is.na(tmp) && tmp > 0 && tmp <= 1) { FDR_THR <- tmp; break }
    cat("Invalid input. Enter a number in (0, 1].\n")
  }

  # Number of permutations
  repeat {
    n_perm_in <- readline(paste0("Permutations per direction (default ", N_PERM, "): "))
    if (nchar(trimws(n_perm_in)) == 0) break
    tmp <- suppressWarnings(as.integer(n_perm_in))
    if (!is.na(tmp) && tmp >= 1) { N_PERM <- tmp; break }
    cat("Invalid input. Please enter a positive integer.\n")
  }

  # Random seed
  seed_in <- readline(paste0("Random seed (default ", RANDOM_SEED, "): "))
  if (nchar(trimws(seed_in)) > 0) {
    tmp <- suppressWarnings(as.integer(seed_in))
    if (!is.na(tmp)) RANDOM_SEED <- tmp
  }
}

cat("FDR threshold:      ", FDR_THR, "\n")
cat("Permutations:       ", N_PERM, "\n")
cat("Random seed:        ", RANDOM_SEED, "\n")

# ============================================================================
# INPUT FILE
# ============================================================================

files <- list.files()

select_file <- function(files, prompt_msg) {
  print(data.frame(index = seq_along(files), filename = files))
  repeat {
    val <- suppressWarnings(as.integer(readline(prompt_msg)))
    if (!is.na(val) && val >= 1 && val <= length(files)) return(files[val])
    cat("Invalid input. Please enter a number between 1 and", length(files), "\n")
  }
}

if (nchar(trimws(RES_FILE)) == 0) {
  RES_FILE <- files[grep("DESeq2_results_test\\.tsv", files, ignore.case = TRUE)][1]
}

if (is.na(RES_FILE) || nchar(trimws(RES_FILE)) == 0) {
  cat("DESeq2_results_test.tsv not found in the working directory.\n")
  if (INTERACTIVE_MODE) {
    RES_FILE <- select_file(files, "Select results file #: ")
  } else {
    stop("RES_FILE not set and DESeq2_results_test.tsv not found in working directory.")
  }
} else if (INTERACTIVE_MODE) {
  if (toupper(readline(paste("Use", RES_FILE, "? (Y/N): "))) != "Y") {
    RES_FILE <- select_file(files, "Select results file #: ")
  }
}

res <- read.delim(RES_FILE, check.names = FALSE)
cat("Rows loaded:", nrow(res), "\n")

need_cols  <- c("GuideID", "Coordinate", "type", "adj-p deplete", "adj-p enrich")
missing_in <- setdiff(need_cols, colnames(res))
if (length(missing_in) > 0) {
  stop("Missing required columns in input TSV: ", paste(missing_in, collapse = ", "))
}

# ============================================================================
# FILTER TO TEST GUIDES
# ============================================================================

res_test <- res[!is.na(res$type) & res$type == "test", ]
if (nrow(res_test) == 0) {
  stop("No rows with type == 'test' found in input. Check the 'type' column.")
}
cat("Test guides retained:", nrow(res_test), "\n")

res_test$Coordinate <- suppressWarnings(as.numeric(as.character(res_test$Coordinate)))
n_before  <- nrow(res_test)
res_test  <- res_test[!is.na(res_test$Coordinate), ]
n_dropped <- n_before - nrow(res_test)
if (n_dropped > 0) cat("Guides dropped (missing coordinate):", n_dropped, "\n")

# Sort guides by coordinate (guide ID breaks ties) so that the permutations,
# and therefore the p-values, do not depend on the row order of the input file.
res_test <- res_test[order(res_test$Coordinate, res_test$GuideID), ]

# ============================================================================
# NEAREST-NEIGHBOR TEST
# ============================================================================

# Median distance from each position to its nearest neighbor.
nn_median <- function(positions) {
  positions <- sort(positions)
  gaps      <- diff(positions)
  median(pmin(c(Inf, gaps), c(gaps, Inf)))
}

run_direction <- function(padj_col) {
  padj      <- suppressWarnings(as.numeric(res_test[[padj_col]]))
  keep      <- !is.na(padj)
  positions <- res_test$Coordinate[keep]
  is_hit    <- padj[keep] < FDR_THR
  n_tested  <- length(positions)
  n_hits    <- sum(is_hit)
  out <- list(n_tested = n_tested, n_hits = n_hits, n_na = sum(!keep),
              observed = NA_real_, p_value = NA_real_,
              null = rep(NA_real_, N_PERM))
  if (n_hits < 2) return(out)

  set.seed(RANDOM_SEED)
  null_stats <- numeric(N_PERM)
  for (i in seq_len(N_PERM)) {
    null_stats[i] <- nn_median(positions[sample.int(n_tested, n_hits)])
  }
  out$observed <- nn_median(positions[is_hit])
  out$p_value  <- (sum(null_stats <= out$observed) + 1) / (N_PERM + 1)
  out$null     <- null_stats
  out
}

directions <- c(depletion = "adj-p deplete", enrichment = "adj-p enrich")
results    <- list()
for (d in names(directions)) {
  cat("\nTesting", d, "hits...\n")
  results[[d]] <- run_direction(directions[[d]])
  r <- results[[d]]
  if (r$n_hits < 2) {
    cat(sprintf("  %d hit(s) at FDR < %g; at least 2 are needed. Skipped.\n",
                r$n_hits, FDR_THR))
  } else {
    cat(sprintf("  Hits: %d of %d guides; observed median NN distance = %.1f bp; p = %.3g\n",
                r$n_hits, r$n_tested, r$observed, r$p_value))
  }
}

# ============================================================================
# OUTPUT
# ============================================================================

dir.create(out_dir, showWarnings = FALSE)
cat("\nOutput folder:", out_dir, "\n")

null_file <- file.path(out_dir, paste0(OUT_PREFIX, "_nn_null_distributions.tsv"))
null_df   <- data.frame(depletion  = results$depletion$null,
                        enrichment = results$enrichment$null)
write.table(null_df, null_file, sep = "\t", quote = FALSE,
            row.names = FALSE, na = "")
cat("Null distributions written:", null_file, "\n")

summary_file <- file.path(out_dir, paste0(OUT_PREFIX, "_nn_summary.tsv"))
summary_df <- data.frame(
  direction                     = names(results),
  hits                          = sapply(results, `[[`, "n_hits"),
  `guides tested`               = sapply(results, `[[`, "n_tested"),
  `observed median NN distance (bp)` = sapply(results, `[[`, "observed"),
  `p-value`                     = sapply(results, `[[`, "p_value"),
  permutations                  = ifelse(sapply(results, `[[`, "n_hits") < 2, NA, N_PERM),
  check.names = FALSE
)
write.table(summary_df, summary_file, sep = "\t", quote = FALSE,
            row.names = FALSE, na = "")
cat("Summary written:", summary_file, "\n")

# ============================================================================
# SESSION & PARAMETER LOG
# ============================================================================

log_file <- file.path(out_dir, paste0(OUT_PREFIX, "_run_parameters.txt"))

result_lines <- unlist(lapply(names(results), function(d) {
  r <- results[[d]]
  if (r$n_hits < 2) {
    sprintf("%-11s hits: %d of %d guides tested — skipped (at least 2 needed)",
            paste0(d, ":"), r$n_hits, r$n_tested)
  } else {
    sprintf("%-11s hits: %d of %d guides tested; observed median NN distance: %.1f bp; p: %.3g",
            paste0(d, ":"), r$n_hits, r$n_tested, r$observed, r$p_value)
  }
}))

log_lines <- c(
  "###############################################################################",
  "# Run parameter log",
  "###############################################################################",
  paste("Timestamp:         ", format(now, "%Y-%m-%d %H:%M:%S")),
  paste("Script version:    ", VERSION),
  paste("Working directory: ", getwd()),
  paste("Input file:        ", RES_FILE),
  paste("Input MD5:         ", unname(tools::md5sum(RES_FILE))),
  paste("Output prefix:     ", OUT_PREFIX),
  paste("Output directory:  ", out_dir),
  "",
  "--- Analysis parameters ---",
  paste("FDR threshold:     ", FDR_THR),
  paste("Permutations:      ", N_PERM),
  paste("Random seed:       ", RANDOM_SEED),
  "",
  "--- Guide counts ---",
  paste("Total rows in input:                ", nrow(res)),
  paste("Test guides (type=='test'):         ", n_before),
  paste("Guides dropped (missing coordinate):", n_dropped),
  paste("Guides dropped (missing adj-p), depletion: ", results$depletion$n_na),
  paste("Guides dropped (missing adj-p), enrichment:", results$enrichment$n_na),
  "",
  "--- Results ---",
  result_lines,
  "",
  "--- R session info ---",
  capture.output(sessionInfo())
)

writeLines(log_lines, log_file)
cat("Parameter log written:", log_file, "\n")

cat("\nDONE\n")
