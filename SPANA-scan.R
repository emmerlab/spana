# SPANA-scan: Saturating Perturbation Analysis of Noncoding Activity
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
# SPANA-scan — sliding window analysis
#
# DESCRIPTION:
#   Aggregates per-guide DESeq2 statistics over sliding genomic windows and
#   computes empirical p-values against a size-matched permutation null.
#   Window means use apeglm-shrunken log2 fold changes, which reduce noise
#   from low-count guides.
#
# REQUIREMENTS:
#   Base R only (the 'parallel' package ships with R).
#
# EXPECTED INPUT (DESeq2_results_test.tsv from SPANA-null):
#   Required columns:
#     GuideID       — unique guide identifier
#     Coordinate    — 1-based genomic coordinate of guide cut site (integer)
#     type          — guide class; only rows where type == "test" are analysed
#     apeglm L2FC   — apeglm-shrunk log2 fold change (experimental / control,
#                     as in the DESeq2 output)
#   Notes:
#     - All coordinates are assumed to be on a single chromosome per run.
#     - Guides with NA coord or NA metric values are dropped.
#     - Control guides (type != "test") are excluded; the permutation null is
#       drawn from test guides only.
#
# NULL MODEL:
#   For each observed window size N, N_SIM window means are simulated by
#   sampling N scores at random (with replacement) from the full test-guide
#   score distribution. Rather than storing all N_SIM simulated means,
#   enrichment and depletion exceedance counts are accumulated incrementally
#   in chunks of SIM_CHUNK_SIZE rows. This keeps peak memory at
#   ~SIM_CHUNK_SIZE * max(N) * 8 bytes regardless of N_SIM, enabling very
#   large simulation counts (e.g. 1e8) without memory constraints.
#   Empirical p-values are the proportion of simulated means >= (enrichment)
#   or <= (depletion) the observed window mean, floored at 1/N_SIM.
#   Reproducibility: guides are sorted by coordinate before simulation, and
#   the simulations for each window size N use the seed RANDOM_SEED + N, so
#   results are identical for a given seed regardless of input row order or
#   the number of cores.
#   NOTE: The window loop is O(n^2) in the number of test guides. For typical
#   screen libraries (~10,000–50,000 guides) runtime is acceptable, but very
#   large libraries may require vectorisation or tiling optimisation.
#
# OUTPUT (apeglm L2FC windows, labelled "LFCshr" in file names; written to a
# timestamped subdirectory):
#   - <prefix>_window_results_LFCshr.tsv     : per-window statistics table
#       (GuideID, Coordinate, window guides, window mean apeglm L2FC,
#        p deplete, adj-p deplete, p enrich, adj-p enrich; BH-adjusted)
#   - WIG files (–log10 p or raw score, variableStep format):
#       <prefix>_window_empirical_pEnrich_LFCshr.wig
#       <prefix>_window_empirical_pDeplete_LFCshr.wig
#       <prefix>_window_mean_LFCshr.wig
#   - BED peak files (overlapping significant windows merged), one per
#     direction and FDR threshold:
#       <prefix>_LFCshr_enrichment_BH_FDR<threshold>_peaks.bed
#       <prefix>_LFCshr_depletion_BH_FDR<threshold>_peaks.bed
#     Each window is written as a WINDOW_BP-wide interval centred on its
#     guide, which is 1 bp wider than the span of cut sites a window can
#     include (guides within ±(WINDOW_BP/2 - 1) bp; e.g. ±49 bp for 100 bp).
#   - <prefix>_run_parameters.txt            : full parameter and session log
#
# WORKING DIRECTORY:
#   Set your working directory to the DESeq2 output folder before sourcing
#   this script, using either the RStudio Files panel (More -> Set As Working
#   Directory) or setwd() in the console. The script will use getwd() if
#   WORK_DIR is left as "".
#
# PARALLELISATION:
#   The null simulation (the main computational bottleneck) is parallelised
#   across unique window sizes N using the 'parallel' package, which ships
#   with base R and requires no additional installation.
#
#   On Linux/macOS (including HPC clusters): mclapply() is used with forking.
#   On Windows: falls back automatically to sequential lapply() since forking
#   is not supported — the script remains fully functional, just single-core.
#
#   On HPC systems, N_CORES is set automatically using nproc, which reports
#   the number of cores allocated to your job rather than the total node
#   cores. This prevents oversubscription on shared nodes. On a computing
#   cluster, request the number of cores you want for the job; the script
#   uses one fewer than the number allocated. On a local machine, the default
#   (detected cores minus 1) is appropriate.
#
#   To bypass all readline() prompts, pre-fill the PARAMETERS block below
#   and set INTERACTIVE_MODE <- FALSE. All prompts will be skipped and the
#   hardcoded values used directly.
###############################################################################


# ============================================================================
# PARAMETERS
# Pre-fill these and set INTERACTIVE_MODE <- FALSE to run non-interactively
# (e.g. on an HPC cluster or within a pipeline).
# ============================================================================

INTERACTIVE_MODE <- TRUE   # Set FALSE to skip all prompts

# --- Parallelisation ---
# Number of cores for null simulation. Default uses all available cores minus
# one to leave the system responsive. On HPC, this is based on the number
# of cores allocated to your session or job; override it here if needed.
# On Windows this setting is ignored and the script runs sequentially.
# On HPC systems, nproc reports allocated cores rather than total node cores.
# Falls back to detectCores() on systems where nproc is unavailable.
nproc_out   <- suppressWarnings(tryCatch(as.integer(system("nproc", intern = TRUE)), error = function(e) NA_integer_))
detected_cores <- suppressWarnings(parallel::detectCores())
if (is.na(detected_cores) || detected_cores < 1L) detected_cores <- 1L  # detection failed: use 1 core
N_CORES     <- if (!is.na(nproc_out) && nproc_out >= 1L) max(1L, nproc_out - 1L) else max(1L, detected_cores - 1L)

# --- Paths ---
WORK_DIR    <- ""          # Working directory; set this before sourcing or use
                           # setwd() in the console. Leave "" to use getwd().
RES_FILE    <- ""          # Input TSV filename (leave "" to auto-detect)
OUT_PREFIX  <- ""          # Output file prefix (e.g. "MyScreen_chr3")

# --- Analysis parameters ---
WINDOW_BP   <- 100L        # Sliding window width in base pairs (must be even)
N_SIM       <- 100000L     # Permutation simulations per window size N
FDR_THRS    <- 0.01        # Peak-calling BH-adjusted p-value threshold
                             # (a BED file is written for each threshold given)
RANDOM_SEED <- 42L         # Random seed for reproducibility

# --- Genome annotation ---
CHROM       <- ""          # Chromosome (e.g. "chr3", "3", "X"); required for the
                           # WIG and BED outputs. Leave "" to be prompted.

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
# OUTPUT DIRECTORY
# ============================================================================

if (INTERACTIVE_MODE && nchar(trimws(OUT_PREFIX)) == 0) {
  OUT_PREFIX <- readline("Output prefix: ")
  if (nchar(trimws(OUT_PREFIX)) == 0) stop("Output prefix cannot be empty.")
}

now    <- Sys.time()
ts     <- format(now, "%Y_%m%d_%H%M")
out_dir <- paste0(ts, "_window_output")

# ============================================================================
# INTERACTIVE OVERRIDES (only when INTERACTIVE_MODE = TRUE)
# ============================================================================

if (INTERACTIVE_MODE) {

  # Number of simulations
  repeat {
    n_sim_input <- readline(paste0("Simulations per N (default ", N_SIM, "): "))
    if (nchar(trimws(n_sim_input)) == 0) break
    tmp <- suppressWarnings(as.integer(n_sim_input))
    if (!is.na(tmp) && tmp >= 1) { N_SIM <- tmp; break }
    cat("Invalid input. Please enter a positive integer.\n")
  }

  # FDR thresholds — prompt accepts comma-separated values (e.g. 0.05,0.01)
  repeat {
    fdr_in <- readline(paste0("Peak-calling FDR thresholds, comma-separated (default ",
                               paste(FDR_THRS, collapse = ","), "): "))
    if (nchar(trimws(fdr_in)) == 0) break
    vals <- suppressWarnings(as.numeric(strsplit(trimws(fdr_in), "[,\\s]+", perl = TRUE)[[1]]))
    if (!any(is.na(vals)) && all(vals > 0) && all(vals <= 1)) { FDR_THRS <- sort(vals); break }
    cat("Invalid input. Enter one or more numbers in (0, 1], comma-separated.\n")
  }

  # Chromosome (required: used to label the WIG tracks and BED peak files)
  while (nchar(trimws(CHROM)) == 0) {
    CHROM <- readline("Chromosome of the screened locus (e.g. 19, X, chrY): ")
    if (nchar(trimws(CHROM)) == 0) cat("A chromosome is required.\n")
  }

  # Random seed
  seed_in <- readline(paste0("Random seed (default ", RANDOM_SEED, "): "))
  if (nchar(trimws(seed_in)) > 0) {
    tmp <- suppressWarnings(as.integer(seed_in))
    if (!is.na(tmp)) RANDOM_SEED <- tmp
  }

}

cat("Simulations per N:  ", N_SIM, "\n")
cat("Cores for sim:      ", N_CORES, "\n")
cat("FDR thresholds:     ", paste(FDR_THRS, collapse = ", "), "\n")
cat("Random seed:        ", RANDOM_SEED, "\n")
cat("Window width (bp):  ", WINDOW_BP, "\n")

# ============================================================================
# DERIVED SETTINGS
# ============================================================================

set.seed(RANDOM_SEED)

# Validate WINDOW_BP — must be a positive even integer for the half_window
# offset to work correctly with integer coordinates.
if (WINDOW_BP %% 2 != 0) stop("WINDOW_BP must be even (current value: ", WINDOW_BP, ")")
if (WINDOW_BP < 2L)      stop("WINDOW_BP must be >= 2 (current value: ", WINDOW_BP, ")")

# Load parallel package (ships with base R; no installation needed)
library(parallel)

# Detect platform for mclapply vs lapply fallback
on_windows <- .Platform$OS.type == "windows"
if (on_windows) {
  cat("Platform: Windows — parallel simulation not available, running sequentially\n")
  N_CORES <- 1L
} else {
  avail <- if (!is.na(nproc_out) && nproc_out >= 1L) nproc_out else detected_cores
  if (N_CORES > avail) {
    cat(sprintf("N_CORES (%d) exceeds available cores (%d) — capping at %d\n",
                N_CORES, avail, avail))
    N_CORES <- avail
  }
  cat(sprintf("Platform: %s — using %d of %d available cores\n",
              Sys.info()[["sysname"]], N_CORES, avail))
}

# half_window: guides within this distance of the focal guide are included.
# For a window of W bp centred on a guide, guides at distance <= (W/2 - 0.5)
# are included (the 0.5 offset ensures the window is symmetric around integer
# coordinates without double-counting boundary positions).
half_window <- WINDOW_BP / 2 - 0.5

# p_floor: prevents empirical p = 0, which would produce infinite –log10(p)
# scores in WIG output. Set to the resolution limit of the simulation.
p_floor <- 1 / N_SIM

# Chromosome formatting (a chromosome is required; WIG and BED files are
# always written)
if (nchar(trimws(CHROM)) == 0) {
  stop("CHROM must be set (e.g. CHROM <- \"19\") when INTERACTIVE_MODE is FALSE.")
}
chrom <- trimws(CHROM)
if (!grepl("^chr", chrom)) chrom <- paste0("chr", chrom)
write_wigs <- TRUE
write_beds <- TRUE

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
  cat("No results TSV auto-detected.\n")
  if (INTERACTIVE_MODE) {
    RES_FILE <- select_file(files, "Select results file #: ")
  } else {
    stop("RES_FILE not set and no *results.tsv found in working directory.")
  }
} else if (INTERACTIVE_MODE) {
  if (toupper(readline(paste("Use", RES_FILE, "? (Y/N): "))) != "Y") {
    RES_FILE <- select_file(files, "Select results file #: ")
  }
}

res <- read.delim(RES_FILE, check.names = FALSE)
cat("Rows loaded:", nrow(res), "\n")

# Map the DESeq2 results-table column names to the names used in this script.
col_map    <- c(guide = "GuideID", coord = "Coordinate", log2FC_shrink = "apeglm L2FC")
missing_in <- setdiff(c(col_map, "type"), colnames(res))
if (length(missing_in) > 0) {
  stop("Missing required columns in input TSV: ", paste(missing_in, collapse = ", "))
}
names(res)[match(col_map, names(res))] <- names(col_map)

# ============================================================================
# FILTER TO TEST GUIDES
# ============================================================================

res_test <- res[!is.na(res$type) & res$type == "test", ]
cat("Test guides retained:", nrow(res_test), "\n")

n_before  <- nrow(res_test)
res_test  <- res_test[!is.na(res_test$coord), ]
n_dropped <- n_before - nrow(res_test)
if (n_dropped > 0) cat("Guides dropped (missing coord):", n_dropped, "\n")

res_test$coord <- as.numeric(as.character(res_test$coord))

# Sort guides by coordinate (guide ID breaks ties) so that the permutation
# null, and therefore the empirical p-values, do not depend on the row order
# of the input file.
res_test <- res_test[order(res_test$coord, res_test$guide), ]

# Required columns check
need_cols <- c("guide", "coord", "log2FC_shrink")
missing   <- setdiff(need_cols, colnames(res_test))
if (length(missing) > 0) {
  stop("Missing required columns in input TSV: ", paste(missing, collapse = ", "))
}

coords <- res_test$coord

# ============================================================================
# INPUT VALIDATION
# ============================================================================

# At least some test guides must remain after filtering
if (nrow(res_test) == 0) {
  stop("No rows with type == 'test' found in input. Check the 'type' column.")
}

# coord must be coercible to finite numeric values
if (all(is.na(res_test$coord))) {
  stop("All 'coord' values are NA after filtering. Check the coordinate column.")
}

# log2FC_shrink must be numeric
if (!is.numeric(res_test$log2FC_shrink) &&
    all(is.na(suppressWarnings(as.numeric(as.character(res_test$log2FC_shrink)))))) {
  stop("'log2FC_shrink' column contains no numeric values. Check DESeq2 output.")
}

# All validation passed — safe to create the output directory now
dir.create(out_dir, showWarnings = FALSE)
cat("Output folder:", out_dir, "\n")

# ============================================================================
# SESSION & PARAMETER LOG
# ============================================================================

log_file <- file.path(out_dir, paste0(OUT_PREFIX, "_run_parameters.txt"))

log_lines <- c(
  "###############################################################################",
  "# Run parameter log",
  "###############################################################################",
  paste("Timestamp:         ", format(now, "%Y-%m-%d %H:%M:%S")),
  paste("Script version:    ", VERSION),
  paste("Working directory: ", getwd()),
  paste("Input file:        ", RES_FILE),
  paste("Output prefix:     ", OUT_PREFIX),
  paste("Output directory:  ", out_dir),
  "",
  "--- Analysis parameters ---",
  paste("Window width (bp): ", WINDOW_BP),
  paste("half_window:       ", half_window),
  paste("N simulations:     ", N_SIM),
  paste("N cores:           ", N_CORES),
  paste("FDR thresholds:    ", paste(FDR_THRS, collapse = ", ")),
  paste("Random seed:       ", RANDOM_SEED),
  paste("p_floor:           ", p_floor),
  paste("Chromosome:        ", chrom),
  "",
  "--- Guide counts ---",
  paste("Total rows in input:          ", nrow(res)),
  paste("Test guides (type=='test'):   ", nrow(res_test)),
  paste("Guides dropped (missing coord):", n_dropped),
  "",
  "--- R session info ---",
  capture.output(sessionInfo())
)

writeLines(log_lines, log_file)
cat("Parameter log written:", log_file, "\n")

# ============================================================================
# HELPERS
# ============================================================================

# Write a variableStep WIG track using –log10(p) as the score.
write_wig_log <- function(df, coord_col, p_col, name) {
  scores <- -log10(pmax(df[[p_col]], 1e-300))
  keep   <- !is.na(df[[coord_col]]) & !is.na(scores) & is.finite(scores)
  if (sum(keep) == 0) {
    cat("No valid rows for track:", name, "— skipping\n")
    return(invisible(NULL))
  }
  file_path <- file.path(out_dir, paste0(OUT_PREFIX, "_", name, ".wig"))
  writeLines(
    c(
      paste0("track type=wiggle_0 name=\"", name, "\""),
      paste0("variableStep chrom=", chrom, " span=1"),
      paste(df[[coord_col]][keep], scores[keep])
    ),
    file_path
  )
  cat("WIG written:", file_path, "(", sum(keep), "rows )\n")
}

# Write a variableStep WIG track using raw values as the score.
write_wig_raw <- function(df, coord_col, val_col, name) {
  vals <- df[[val_col]]
  keep <- !is.na(df[[coord_col]]) & !is.na(vals) & is.finite(vals)
  if (sum(keep) == 0) {
    cat("No valid rows for track:", name, "— skipping\n")
    return(invisible(NULL))
  }
  file_path <- file.path(out_dir, paste0(OUT_PREFIX, "_", name, ".wig"))
  writeLines(
    c(
      paste0("track type=wiggle_0 name=\"", name, "\""),
      paste0("variableStep chrom=", chrom, " span=1"),
      paste(df[[coord_col]][keep], vals[keep])
    ),
    file_path
  )
  cat("WIG written:", file_path, "(", sum(keep), "rows )\n")
}

# Build a BED-format data.frame of significant windows.
# Window boundaries are derived from half_window so they stay consistent
# if WINDOW_BP is changed above.
make_windows_bed <- function(df, padj_col, mean_col, fdr_thr = NULL) {
  if (is.null(fdr_thr)) fdr_thr <- min(FDR_THRS)
  keep <- !is.na(df$coord) &
          !is.na(df[[padj_col]]) &
          is.finite(df[[padj_col]]) &
          df[[padj_col]] <= fdr_thr
  sig <- df[keep, ]
  if (nrow(sig) == 0) return(NULL)

  # Convert 1-based window centre to 0-based half-open BED coordinates.
  # Window spans [coord - half_window - 0.5, coord + half_window + 0.5)
  # which rounds to [coord - WINDOW_BP/2, coord + WINDOW_BP/2). This interval
  # is 1 bp wider than the span of cut sites a window can include (see the
  # OUTPUT section of the header).
  half_bp   <- as.integer(WINDOW_BP / 2L)
  start_0   <- pmax(as.integer(round(sig$coord)) - half_bp, 0L)
  end_0     <- as.integer(round(sig$coord)) + half_bp

  data.frame(
    chrom     = chrom,
    start     = start_0,
    end       = end_0,
    coord     = sig$coord,
    meanScore = sig[[mean_col]],
    padj      = sig[[padj_col]],
    stringsAsFactors = FALSE
  )
}

# Merge overlapping BED windows and summarise each merged peak.
merge_intervals_with_summary <- function(bed_df) {
  bed_df    <- bed_df[order(bed_df$start, bed_df$end), ]
  peaks     <- list()
  cur_start <- bed_df$start[1]
  cur_end   <- bed_df$end[1]
  idxs      <- 1L

  flush_peak <- function(idxs, cur_start, cur_end) {
    sc <- bed_df$meanScore[idxs]
    pj <- bed_df$padj[idxs]
    data.frame(
      chrom            = bed_df$chrom[1],
      start            = cur_start,
      end              = cur_end,
      meanScore        = mean(sc, na.rm = TRUE),
      meanNegLog10Padj = mean(-log10(pmax(pj, 1e-300)), na.rm = TRUE),
      n_windows        = length(idxs)
    )
  }

  # Guard: a single significant window needs no merging — return it directly.
  # Without this, `for (i in 2:1)` would iterate backwards (R does not treat
  # this as a no-op), corrupting the output.
  if (nrow(bed_df) == 1L) {
    outp <- flush_peak(1L, cur_start, cur_end)
    rownames(outp) <- NULL
    return(outp)
  }

  for (i in 2:nrow(bed_df)) {
    if (bed_df$start[i] <= cur_end) {
      cur_end <- max(cur_end, bed_df$end[i])
      idxs    <- c(idxs, i)
    } else {
      peaks[[length(peaks) + 1L]] <- flush_peak(idxs, cur_start, cur_end)
      cur_start <- bed_df$start[i]
      cur_end   <- bed_df$end[i]
      idxs      <- i
    }
  }
  peaks[[length(peaks) + 1L]] <- flush_peak(idxs, cur_start, cur_end)

  outp <- do.call(rbind, peaks)
  rownames(outp) <- NULL
  outp
}

# Write a merged-peaks BED file (BED6 + extra annotation columns).
write_bed_peaks <- function(peaks_df, name) {
  if (is.null(peaks_df) || nrow(peaks_df) == 0) {
    cat("No peaks to write for:", name, "\n")
    return(invisible(NULL))
  }
  # BED score field: scale meanNegLog10Padj to 0–1000 (UCSC convention)
  score     <- pmin(as.integer(round(peaks_df$meanNegLog10Padj * 100)), 1000L)
  peak_name <- paste0(name, "_", seq_len(nrow(peaks_df)))

  bed_out <- data.frame(
    chrom            = peaks_df$chrom,
    start            = peaks_df$start,
    end              = peaks_df$end,
    name             = peak_name,
    score            = score,
    strand           = ".",
    meanScore        = peaks_df$meanScore,
    meanNegLog10Padj = peaks_df$meanNegLog10Padj,
    n_windows        = peaks_df$n_windows
  )

  bed_file <- file.path(out_dir, paste0(OUT_PREFIX, "_", name, "_peaks.bed"))
  write.table(bed_out, bed_file, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  cat("BED written:", bed_file, "(", nrow(bed_out), "peaks )\n")
}

# ============================================================================
# CORE WINDOW ANALYSIS FUNCTION (called once per metric)
# ============================================================================

run_window_analysis <- function(metric_label, score_vec, mean_col_label) {

  # --- Filter to guides with a valid score for this metric ---
  keep_metric <- !is.na(score_vec) & is.finite(score_vec) & !is.na(coords)
  dfm      <- res_test[keep_metric, ]
  coords_m <- dfm$coord
  scores_m <- score_vec[keep_metric]

  cat("\n------------------------------------------\n")
  cat("Window analysis:", metric_label, "\n")
  cat("Guides with valid scores:", nrow(dfm), "\n")

  # --- Sliding window summaries (O(n^2); see header note on large libraries) ---
  n_guides   <- integer(nrow(dfm))
  mean_score <- numeric(nrow(dfm))

  for (i in seq_len(nrow(dfm))) {
    in_window     <- abs(coords_m - coords_m[i]) <= half_window
    n_guides[i]   <- sum(in_window)
    mean_score[i] <- mean(scores_m[in_window])
  }

  cat("Unique window sizes (N):", length(unique(n_guides)), "\n")
  cat("Window N range: [", min(n_guides), ",", max(n_guides), "]\n")
  cat("------------------------------------------\n")

  # --- Incremental permutation null and empirical p-values ---
  #
  # For each unique window size N, N_SIM window means are simulated by
  # sampling N scores at random (with replacement) from the observed
  # score distribution. Rather than storing all N_SIM means, enrichment
  # and depletion exceedance counts are accumulated incrementally in
  # chunks of SIM_CHUNK_SIZE rows. This keeps peak memory at
  # ~SIM_CHUNK_SIZE * max(N) * 8 bytes regardless of N_SIM, and null
  # storage at O(1) per window size rather than O(N_SIM).
  #
  # Results are mathematically identical to the full-matrix approach:
  #   p_enrich = (simulated means >= observed) / N_SIM
  #   p_deplete = (simulated means <= observed) / N_SIM
  # Both are floored at p_floor = 1/N_SIM (simulation resolution limit).
  #
  # Parallelisation: each unique N is independent. On Linux/macOS,
  # mclapply() forks N_CORES workers. On Windows, falls back to lapply().
  unique_n       <- sort(unique(n_guides))
  SIM_CHUNK_SIZE <- 10000L

  cat(sprintf("  Simulating %d unique window sizes using %d core(s)...\n",
              length(unique_n), N_CORES))

  # For each unique N, return a two-element vector:
  #   [1] number of simulated means >= each observed window mean (summed)
  #   [2] number of simulated means <= each observed window mean (summed)
  # obs_by_n: named list mapping N -> vector of observed means for windows
  #           of that size (needed inside the worker)
  obs_by_n <- lapply(
    setNames(unique_n, as.character(unique_n)),
    function(n) mean_score[n_guides == n]
  )

  sim_one_n <- function(n, obs_by_n, scores_m, N_SIM, SIM_CHUNK_SIZE, p_floor, seed) {
    # Each window size sets its own seed (run seed + window size), so the
    # random draws, and therefore the p-values, are identical regardless of
    # the number of cores or the order in which window sizes are processed.
    set.seed(seed + n)
    obs_vec    <- obs_by_n[[as.character(n)]]
    n_obs      <- length(obs_vec)
    en_counts  <- integer(n_obs)   # counts of sim means >= obs
    de_counts  <- integer(n_obs)   # counts of sim means <= obs
    n_chunks   <- ceiling(N_SIM / SIM_CHUNK_SIZE)
    idx_start  <- 1L
    # Report progress every 10% of chunks
    report_at  <- unique(round(seq(0.1, 1.0, by = 0.1) * n_chunks))
    for (chunk in seq_len(n_chunks)) {
      chunk_rows  <- min(SIM_CHUNK_SIZE, N_SIM - idx_start + 1L)
      chunk_means <- rowMeans(matrix(
        sample(scores_m, size = chunk_rows * n, replace = TRUE),
        nrow = chunk_rows,
        ncol = n
      ))
      # For each observed window of this size, count exceedances
      for (j in seq_len(n_obs)) {
        en_counts[j] <- en_counts[j] + sum(chunk_means >= obs_vec[j])
        de_counts[j] <- de_counts[j] + sum(chunk_means <= obs_vec[j])
      }
      idx_start <- idx_start + chunk_rows
      if (chunk %in% report_at) {
        message(sprintf("  N=%d: %d%% complete (%d/%d chunks)",
                        n, as.integer(round(100 * chunk / n_chunks)),
                        chunk, n_chunks))
      }
    }
    list(en = en_counts, de = de_counts)
  }

  # Pass all required variables explicitly as arguments to avoid lexical
  # scoping failures in forked mclapply workers.
  sim_args <- list(
    obs_by_n       = obs_by_n,
    scores_m       = scores_m,
    N_SIM          = N_SIM,
    SIM_CHUNK_SIZE = SIM_CHUNK_SIZE,
    p_floor        = p_floor,
    seed           = RANDOM_SEED
  )

  result_list <- if (on_windows || N_CORES == 1L) {
    lapply(unique_n, function(n) do.call(sim_one_n, c(list(n = n), sim_args)))
  } else {
    parallel::mclapply(unique_n,
                       function(n) do.call(sim_one_n, c(list(n = n), sim_args)),
                       mc.cores = N_CORES)
  }
  names(result_list) <- as.character(unique_n)

  cat("  Simulation complete\n")

  # --- Assemble empirical p-values from exceedance counts ---
  win_p_en <- numeric(nrow(dfm))
  win_p_de <- numeric(nrow(dfm))

  # For each unique N, map results back to the correct window indices
  for (n_char in names(result_list)) {
    idx_in_dfm  <- which(n_guides == as.integer(n_char))
    res         <- result_list[[n_char]]
    win_p_en[idx_in_dfm] <- pmax(res$en / N_SIM, p_floor)
    win_p_de[idx_in_dfm] <- pmax(res$de / N_SIM, p_floor)
  }

  # --- Multiple testing correction (Benjamini-Hochberg) ---
  non_na     <- !is.na(win_p_en)
  padj_en_BH <- rep(NA_real_, length(win_p_en))
  padj_de_BH <- rep(NA_real_, length(win_p_de))

  padj_en_BH[non_na] <- p.adjust(win_p_en[non_na], method = "BH")
  padj_de_BH[non_na] <- p.adjust(win_p_de[non_na], method = "BH")

  # --- Assemble output table ---
  outm <- data.frame(
    guide          = dfm$guide,
    coord          = coords_m,
    window_n       = n_guides,
    window_mean    = mean_score,
    p_enrich       = win_p_en,
    padj_enrich_BH = padj_en_BH,
    p_deplete      = win_p_de,
    padj_deplete_BH = padj_de_BH,
    stringsAsFactors = FALSE
  )

  outm <- outm[order(outm$coord), ]
  colnames(outm)[colnames(outm) == "window_mean"] <- mean_col_label

  # --- Summary (all FDR thresholds) ---
  cat("\nSUMMARY:", metric_label, "\n")
  cat("  Windows computed:                   ", nrow(outm), "\n")
  cat("  Mean guides per window:             ", round(mean(outm$window_n), 1), "\n")
  cat("  Mean window score:                  ", round(mean(outm[[mean_col_label]]), 5), "\n")
  for (thr in FDR_THRS) {
    cat(sprintf("  Enriched  (BH < %.2f):              %d\n", thr,
                sum(outm$padj_enrich_BH  < thr, na.rm = TRUE)))
    cat(sprintf("  Depleted  (BH < %.2f):              %d\n", thr,
                sum(outm$padj_deplete_BH < thr, na.rm = TRUE)))
  }

  # --- Write window results TSV ---
  # Column names and order follow the supplementary window tables; Coordinate
  # is the cut site coordinate of the guide at the center of each window.
  out_tab <- data.frame(
    GuideID       = outm$guide,
    Coordinate    = outm$coord,
    window_guides = outm$window_n,
    window_mean   = outm[[mean_col_label]],
    p_deplete     = outm$p_deplete,
    adjp_deplete  = outm$padj_deplete_BH,
    p_enrich      = outm$p_enrich,
    adjp_enrich   = outm$padj_enrich_BH,
    stringsAsFactors = FALSE
  )
  colnames(out_tab) <- c("GuideID", "Coordinate", "window guides",
                         "window mean apeglm L2FC", "p deplete", "adj-p deplete",
                         "p enrich", "adj-p enrich")
  out_file <- file.path(
    out_dir, paste0(OUT_PREFIX, "_window_results_", metric_label, ".tsv")
  )
  write.table(out_tab, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Window results written:", out_file, "\n")

  # --- WIG tracks (unadjusted empirical p-values and raw window mean) ---
  if (write_wigs) {
    write_wig_log(outm, "coord", "p_enrich",
                  paste0("window_empirical_pEnrich_",  metric_label))
    write_wig_log(outm, "coord", "p_deplete",
                  paste0("window_empirical_pDeplete_", metric_label))
    write_wig_raw(outm, "coord", mean_col_label,
                  paste0("window_mean_",               metric_label))
  }

  # --- BED peak files (all FDR thresholds, BH; enrichment and depletion) ---
  if (write_beds) {
    for (thr in FDR_THRS) {
      thr_label <- gsub("\\.", "", sprintf("%.2f", thr))  # e.g. "0.05" -> "005"
      for (m in c("BH")) {
        bed_en <- make_windows_bed(outm, paste0("padj_enrich_",  m), mean_col_label,
                                   fdr_thr = thr)
        if (!is.null(bed_en)) {
          write_bed_peaks(merge_intervals_with_summary(bed_en),
                          paste0(metric_label, "_enrichment_", m, "_FDR", thr_label))
        } else {
          cat(sprintf("  No enrichment windows pass FDR %.2f for %s %s\n",
                      thr, metric_label, m))
        }

        bed_de <- make_windows_bed(outm, paste0("padj_deplete_", m), mean_col_label,
                                   fdr_thr = thr)
        if (!is.null(bed_de)) {
          write_bed_peaks(merge_intervals_with_summary(bed_de),
                          paste0(metric_label, "_depletion_", m, "_FDR", thr_label))
        } else {
          cat(sprintf("  No depletion windows pass FDR %.2f for %s %s\n",
                      thr, metric_label, m))
        }
      }
    }
  }

  invisible(outm)
}

# ============================================================================
# RUN WINDOW ANALYSIS
# ============================================================================

# Shrunken log2FC windows (apeglm posterior median), which reduce noise from
# low-count guides.
run_window_analysis(
  metric_label   = "LFCshr",
  score_vec      = as.numeric(as.character(res_test$log2FC_shrink)),
  mean_col_label = "window_mean_log2FC_shrink"
)

cat("\nDONE\n")
