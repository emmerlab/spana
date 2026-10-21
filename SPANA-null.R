# SPANA-null: Saturating Perturbation Analysis of Noncoding Activity
#
# Copyright (C) 2026 The Regents of the University of Michigan
# Emmer Lab
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 2 of the License, or (at your option)
# any later version. It is distributed WITHOUT ANY WARRANTY; see the LICENSE
# file for details.

VERSION <- "1.0.0"  # script version; recorded in the run log

# =============================================================================
# SPANA-null — null-set normalization and guide-level DESeq2 analysis
# =============================================================================
#
# PURPOSE
#   Differential abundance analysis of individual guide RNAs from a pooled
#   CRISPR screen with two FACS-sorted populations (experimental vs control).
#   Designed for noncoding screens in which asymmetric depletion or
#   enrichment of guides targeting functional elements can violate the
#   assumptions of standard median-of-ratios normalization.
#
# NORMALIZATION
#   Size factors are estimated from an empirically defined null guide set
#   using iterative sigma clipping of the log2 fold change (LFC) distribution.
#   Guides beyond ±k SD of the mean are excluded, the mean and SD are
#   recomputed on the remaining guides, and this is repeated until membership
#   stabilizes. Size factors are then reestimated, LFC values recomputed, and
#   the clipping procedure repeated until null set membership is stable across
#   passes (Jaccard similarity >= 0.98). This approach does not require
#   pre-specifying any genomic region or control guides as a normalization
#   reference.
#
# SIGNIFICANCE TESTING
#   One-sided p-values for depletion and enrichment are derived independently
#   from the DESeq2 Wald statistic (standard normal under the null).
#   Benjamini-Hochberg (BH) adjusted p-values are reported for each direction.
#   Shrunken LFC estimates are computed using apeglm.
#
# INPUTS (prompted interactively at runtime)
#   - Counts file:       tab-delimited, guides as rows, samples as columns.
#                        The first column holds guide IDs; every other column
#                        is read as a sample, so the file should contain no
#                        other columns.
#   - Annotations file:  tab-delimited, required columns (names matched
#                        without regard to capitalization):
#                          Guide_ID    — must match row names of counts file
#                          Guide_Type  — "test", "pos", or "neg"
#                          Coordinate  — genomic cut site coordinate (integer);
#                                        may be blank for guides without one
#                        optional columns carried through to output as-is.
#                        Numbers written with thousands separators (e.g.
#                        "11,039,456") are read as plain numbers.
#   - Sample assignment: detected automatically from column names or entered
#                        manually; user confirms control vs experimental bins
#                        and replicate labels
#
# KEY PARAMETERS (all have defaults; k [k_clip in the script], n_passes,
# jaccard_thresh, and min_null_n are set under "Normalization settings" in
# section 5; min_count and min_samples are set in section 4). All are
# recorded in the run log.
#   k              = 2.0   Sigma clipping threshold (SD units).
#   n_passes       = 5     Maximum outer normalization iterations. Convergence
#                          typically achieved in 2-3 passes.
#   jaccard_thresh = 0.98  Convergence criterion for null set membership.
#   min_null_n     = 100   Minimum null set size. Increase k if this is hit.
#   min_count      = 10    Minimum count per guide per sample for pre-filtering.
#   min_samples    = 3     Minimum number of samples meeting min_count.
#
# OUTPUTS (written to a timestamped folder in the working directory)
#   DESeq2_results_test.tsv          — results for test guides (input to the
#                                      window analysis)
#   DESeq2_results_controls.tsv      — results for pos and neg control guides
#   QC_PCA.tsv                       — PCA of VST-normalized counts: PC scores
#                                      per sample, with percent variance
#                                      explained in each column header
#   run_log.txt                      — provenance, parameters, guide accounting,
#                                      normalization passes, results summary
#
# OUTPUT COLUMN DESCRIPTIONS (both results tables)
#   GuideID            guide ID (from counts file row names)
#   type               guide type (test / pos / neg)
#   Coordinate         genomic cut site coordinate (from the annotations file)
#   [annotation_cols_front]  annotation columns, as named in the settings below
#   Raw L2FC           MLE LFC (exp / control), DESeq2 Wald test
#   Raw L2FC SE        standard error of Raw L2FC
#   stat               Wald statistic (standard normal under null)
#   apeglm L2FC        apeglm posterior median LFC (preferred for ranking)
#   apeglm L2FC SE     posterior SD of apeglm L2FC
#   p deplete          one-sided p-value, depletion (left tail of Wald stat)
#   adj-p deplete      BH-adjusted p deplete
#   p enrich           one-sided p-value, enrichment (right tail of Wald stat)
#   adj-p enrich       BH-adjusted p enrich
#   baseMean           mean normalized count across all samples
#   norm-cnt <sample>  normalized count per sample (raw count / size factor);
#                      <sample> is the sample name as written in the counts file
#   L2FC <rep>         per-replicate LFC: log2((norm_exp+1)/(norm_ctrl+1))
#   [annotation_cols_back]   annotation columns, as named in the settings below
#   DESeq2 Null Group  TRUE if guide was in the converged sigma-clipping null
#                      set (test table only)
#   Control guide adjusted p-values are corrected within each control type and
#   are not comparable to test guide adjusted p-values.
#
# DEPENDENCIES
#   install.packages(c("tidyverse", "BiocManager"))
#   BiocManager::install(c("DESeq2", "apeglm"))
#   matrixStats is also used (PCA step); it is installed automatically
#   as a dependency of DESeq2.
# =============================================================================

library(DESeq2)
library(apeglm)
library(tidyverse)

# Set working directory before sourcing, or set here:
# setwd("/path/to/your/directory")

# --- Output table settings ---------------------------------------------------
# Optional annotation columns to include in the results tables, placed by name:
annotation_cols_front <- c("group")               # placed after 'Coordinate'
annotation_cols_back  <- c("LDLR CDS proximity")  # placed after the 'L2FC <rep>' columns

# =============================================================================
# 0. OUTPUT DIRECTORY AND RUN LOG
# =============================================================================

timestamp  <- format(Sys.time(), "%Y_%m_%d_%H%M")
output_dir <- paste0(timestamp, "_DESeq2_output")
dir.create(output_dir)
cat("Output directory:", output_dir, "\n\n")

# --- Select input files from working directory -------------------------------

available_files <- list.files(getwd())

cat("Available files in working directory:\n")
for (i in seq_along(available_files)) {
  cat(sprintf("  [%d] %s\n", i, available_files[i]))
}

counts_idx       <- as.integer(readline(prompt = "\nEnter number for counts file: "))
counts_file      <- available_files[counts_idx]

annotations_idx  <- as.integer(readline(prompt = "Enter number for annotations file: "))
annotations_file <- available_files[annotations_idx]

cat(sprintf("\nCounts file:      %s\n", counts_file))
cat(sprintf("Annotations file: %s\n\n", annotations_file))

# --- Attempt to record script path for provenance ---------------------------

script_path <- tryCatch(
  # When the script is run with source(), sys.frame(1)$ofile gives its path.
  # When it is run from the command line with Rscript, that fails and the
  # --file= argument is used instead. If neither is available (e.g. the code
  # is pasted into the console), the path is recorded as unknown.
  normalizePath(sys.frame(1)$ofile),
  error = function(e) {
    tryCatch(
      normalizePath(commandArgs(trailingOnly = FALSE) %>%
                      .[grepl("--file=", .)] %>%
                      gsub("--file=", "", .)),
      error = function(e) character(0)
    )
  }
)
if (length(script_path) == 0 || !nzchar(script_path[1]))
  script_path <- "Unknown (not run with source() or Rscript)"

# --- Write run log -----------------------------------------------------------

log_lines <- c(
  "=============================================================================",
  "RUN LOG",
  "=============================================================================",
  sprintf("Timestamp:          %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("Script:             %s", script_path),
  sprintf("Script version:     %s", VERSION),
  sprintf("Working directory:  %s", normalizePath(getwd())),
  "",
  "Input files:",
  sprintf("  Counts:           %s", normalizePath(counts_file)),
  sprintf("  Annotations:      %s", normalizePath(annotations_file)),
  sprintf("  Counts MD5:       %s", tools::md5sum(counts_file)),
  sprintf("  Annotations MD5:  %s", tools::md5sum(annotations_file)),
  "",
  "R environment:",
  sprintf("  R version:        %s", R.version$version.string),
  sprintf("  DESeq2:           %s", as.character(packageVersion("DESeq2"))),
  sprintf("  apeglm:           %s", as.character(packageVersion("apeglm"))),
  sprintf("  tidyverse:        %s", as.character(packageVersion("tidyverse"))),
  "",
  sprintf("Output directory:   %s", normalizePath(output_dir)),
  "============================================================================="
)

log_path <- file.path(output_dir, "run_log.txt")
writeLines(log_lines, log_path)
cat(paste(log_lines, collapse = "\n"), "\n\n")

# =============================================================================
# 1. LOAD DATA
# =============================================================================

# check.names = FALSE keeps sample and annotation column names exactly as
# written in the input files (spaces and hyphens are not converted to dots),
# so they can be carried through to the output column headers.
counts_raw  <- read.delim(counts_file, row.names = 1, check.names = FALSE)
annotations <- read.delim(annotations_file, check.names = FALSE)

# --- Standardize annotation column names and number formatting --------------
# Required columns are matched without regard to capitalization, so a header
# of "guide_ID" or "GUIDE_TYPE" is read as Guide_ID / Guide_Type.
for (req in c("Guide_ID", "Guide_Type", "Coordinate")) {
  hit <- which(tolower(trimws(colnames(annotations))) == tolower(req))
  if (length(hit) == 0)
    stop(sprintf("Annotations file is missing required column '%s' (any capitalization).", req))
  if (length(hit) > 1)
    stop(sprintf("Annotations file has more than one column matching '%s': %s",
                 req, paste(colnames(annotations)[hit], collapse = ", ")))
  colnames(annotations)[hit] <- req
}

# Numbers saved with thousands separators (e.g. "11,039,456", as exported
# from Excel) are read as text; convert any column whose non-blank values
# are all numbers once commas are removed. Other columns are left unchanged.
strip_commas_if_numeric <- function(x) {
  if (!is.character(x)) return(x)
  y <- gsub(",", "", trimws(x))
  y[y == ""] <- NA
  num <- suppressWarnings(as.numeric(y))
  if (all(is.na(y) | !is.na(num))) num else x
}
annotations <- annotations %>% mutate(across(-Guide_ID, strip_commas_if_numeric))

if (!is.numeric(annotations$Coordinate))
  stop("The Coordinate column contains values that are not numbers. ",
       "Check the annotations file.")

# --- Check Guide_Type labels -------------------------------------------------
# Labels are read without regard to capitalization or surrounding spaces
# ("Test" or "test " is read as "test"). Any other label stops the script,
# because guides with an unrecognized type would otherwise be silently
# excluded from both the analysis and the output.
guide_types <- c("test", "pos", "neg")
annotations$Guide_Type <- tolower(trimws(as.character(annotations$Guide_Type)))
bad_type <- is.na(annotations$Guide_Type) | !annotations$Guide_Type %in% guide_types
if (any(bad_type)) {
  bad_tab <- table(ifelse(is.na(annotations$Guide_Type[bad_type]), "(blank)",
                          annotations$Guide_Type[bad_type]))
  stop(sprintf(
    "Unrecognized Guide_Type value(s) in the annotations file (expected test, pos, or neg):%s",
    paste0("\n  '", names(bad_tab), "': ", as.integer(bad_tab), " guide(s)", collapse = "")
  ))
}

# --- Drop annotated guides that are absent from the counts file -------------
# Guides can be missing from the counts file (e.g. if no reads were counted
# for them). They are removed here, before any analysis, with a warning, and
# the number removed per guide type is recorded in the run log (section 4).
n_annotated_by_type <- table(factor(annotations$Guide_Type, levels = guide_types))
missing_counts      <- setdiff(annotations$Guide_ID, rownames(counts_raw))
n_missing_by_type   <- table(factor(
  annotations$Guide_Type[annotations$Guide_ID %in% missing_counts], levels = guide_types))
if (length(missing_counts) > 0) {
  warning(sprintf(
    "%d annotated guide(s) are absent from the counts file and will be excluded:%s%s",
    length(missing_counts),
    paste0("\n  ", head(missing_counts, 10), collapse = ""),
    if (length(missing_counts) > 10)
      sprintf("\n  ... and %d more", length(missing_counts) - 10) else ""
  ))
  annotations <- annotations[!annotations$Guide_ID %in% missing_counts, , drop = FALSE]
}

# Drop guides present in the counts file but absent from the annotations file,
# with a warning listing the excluded IDs. This handles spike-in guides or
# guides from a different library that were not annotated for this screen.
missing_ann <- setdiff(rownames(counts_raw), annotations$Guide_ID)
if (length(missing_ann) > 0) {
  warning(sprintf(
    "%d guide(s) in counts file have no annotation and will be excluded:%s%s",
    length(missing_ann),
    paste0("\n  ", head(missing_ann, 10), collapse = ""),
    if (length(missing_ann) > 10)
      sprintf("\n  ... and %d more", length(missing_ann) - 10) else ""
  ))
  counts_raw <- counts_raw[rownames(counts_raw) %in% annotations$Guide_ID, , drop = FALSE]
}
cat(sprintf("Guides in counts file: %d total, %d annotated, %d dropped\n\n",
            length(missing_ann) + nrow(counts_raw),
            nrow(counts_raw),
            length(missing_ann)))

# =============================================================================
# 2. FILTER TO TEST GUIDES
#
# Three guide types are expected in the annotations file:
#   "test" — guides targeting the locus of interest; primary analysis set
#   "pos"  — positive controls; not used for normalization; analyzed
#            separately for validation and quality control
#   "neg"  — negative controls; not used for normalization; analyzed
#            separately for validation and quality control
#
# pos and neg guides are processed separately in section 11 using the
# size factors derived from test guide normalization, placing their LFC
# values on the same scale for validation purposes.
# =============================================================================

test_guides <- annotations %>%
  filter(Guide_Type == "test") %>%
  pull(Guide_ID)

counts_test <- counts_raw[test_guides, ]

cat("Test guides:", nrow(counts_test), "\n")
cat("Control guides (pos + neg, analyzed separately):",
    nrow(counts_raw) - nrow(counts_test), "\n\n")

# =============================================================================
# 3. SAMPLE ASSIGNMENT
#
# Sample columns are read from the counts file and assigned to two bins
# (control and experimental) with replicate labels. Assignment is attempted
# automatically by scanning for contrasting keyword pairs in column names;
# the user confirms or overrides the guess. The user then specifies which
# bin is the reference (control). Positive LFC = enriched in experimental.
# =============================================================================

samples <- colnames(counts_test)
cat("Sample columns detected:\n")
for (i in seq_along(samples)) cat(sprintf("  [%d] %s\n", i, samples[i]))
cat("\n")

# --- Helper: infer replicate number from sample name -----------------------
# Tries patterns in order: RepN, RN, _N, .N, trailing integer.
# Returns NULL if no consistent match is found.
guess_replicates <- function(nms) {
  patterns <- c("[Rr]ep?([0-9]+)", "[Rr]([0-9]+)$",
                "_([0-9]+)$",      "\\.([0-9]+)$",
                "([0-9]+)$")
  for (pat in patterns) {
    m <- regmatches(nms, regexpr(pat, nms, perl = TRUE))
    if (length(m) == length(nms) && !any(is.na(m))) {
      nums <- gsub("[^0-9]", "", m)
      if (length(unique(nums)) == length(nms)) return(paste0("Rep", nums))
    }
  }
  return(NULL)
}

# --- Helper: manual bin assignment -----------------------------------------
manual_assign_bins <- function(samples) {
  cat("Manual bin assignment.\n")
  cat("Enter space- or comma-separated indices for the CONTROL bin:\n")
  ctrl_input <- readline(prompt = "  Control indices: ")
  ctrl_idx   <- as.integer(strsplit(trimws(ctrl_input), "[,\\s]+", perl = TRUE)[[1]])

  cat("Enter space- or comma-separated indices for the EXPERIMENTAL bin:\n")
  exp_input <- readline(prompt = "  Experimental indices: ")
  exp_idx   <- as.integer(strsplit(trimws(exp_input), "[,\\s]+", perl = TRUE)[[1]])

  if (any(is.na(ctrl_idx)) || any(is.na(exp_idx)))
    stop("Invalid indices entered during manual bin assignment.")
  if (length(intersect(ctrl_idx, exp_idx)) > 0)
    stop("The same sample cannot be assigned to both bins.")
  if (!all(c(ctrl_idx, exp_idx) %in% seq_along(samples)))
    stop("One or more indices are out of range.")

  list(control = samples[ctrl_idx], exp = samples[exp_idx])
}

# --- Helper: manual replicate label assignment -----------------------------
manual_assign_reps <- function(bin_samples, bin_label) {
  cat(sprintf("Enter replicate labels for %s bin (space-separated,\n", bin_label))
  cat(sprintf("in the same order as: %s):\n", paste(bin_samples, collapse = ", ")))
  rep_input  <- readline(prompt = "  Replicate labels: ")
  rep_labels <- strsplit(trimws(rep_input), "\\s+")[[1]]
  if (length(rep_labels) != length(bin_samples))
    stop(sprintf(
      "Number of labels (%d) does not match number of %s samples (%d).",
      length(rep_labels), bin_label, length(bin_samples)
    ))
  rep_labels
}

# --- Automatic bin detection -----------------------------------------------
# Scans for contrasting keyword pairs common in FACS-sorted screens.
bin_keyword_pairs <- list(
  c("high", "low"), c("hi", "lo"), c("enriched", "depleted"),
  c("top", "bottom"), c("plus", "minus"),
  c("bin1", "bin2"), c("sort1", "sort2"), c("gate1", "gate2")
)

guess_bins <- NULL
for (pair in bin_keyword_pairs) {
  idx_a <- grepl(pair[1], samples, ignore.case = TRUE)
  idx_b <- grepl(pair[2], samples, ignore.case = TRUE)
  if (any(idx_a) && any(idx_b) && !any(idx_a & idx_b)) {
    guess_bins <- list(a = samples[idx_a], b = samples[idx_b],
                       kw_a = pair[1], kw_b = pair[2])
    break
  }
}

# --- Confirm or override bin assignment ------------------------------------
if (!is.null(guess_bins)) {
  cat(sprintf("Detected keyword pair: '%s' / '%s'\n",
              guess_bins$kw_a, guess_bins$kw_b))
  cat(sprintf("  Group A ('%s'): %s\n", guess_bins$kw_a,
              paste(guess_bins$a, collapse = ", ")))
  cat(sprintf("  Group B ('%s'): %s\n\n", guess_bins$kw_b,
              paste(guess_bins$b, collapse = ", ")))
  if (toupper(trimws(readline("Is this bin assignment correct? (Y/N): "))) == "Y") {
    bins_ab <- list(a = guess_bins$a, b = guess_bins$b)
  } else {
    bins_ab <- manual_assign_bins(samples)
    bins_ab <- list(a = bins_ab$control, b = bins_ab$exp)
  }
} else {
  cat("Could not automatically detect bin groupings from sample names.\n")
  bins_ab <- manual_assign_bins(samples)
  bins_ab <- list(a = bins_ab$control, b = bins_ab$exp)
}

# --- Designate control vs experimental -------------------------------------
cat("\nWhich group is the CONTROL (reference) bin?\n")
cat(sprintf("  [1] %s\n", paste(bins_ab$a, collapse = ", ")))
cat(sprintf("  [2] %s\n", paste(bins_ab$b, collapse = ", ")))
ref_choice <- as.integer(trimws(readline(prompt = "Enter 1 or 2: ")))
if (!ref_choice %in% c(1L, 2L)) stop("Invalid choice: enter 1 or 2.")

ctrl_samples <- if (ref_choice == 1L) bins_ab$a else bins_ab$b
exp_samples  <- if (ref_choice == 1L) bins_ab$b else bins_ab$a

cat(sprintf("\nControl bin:      %s\n", paste(ctrl_samples, collapse = ", ")))
cat(sprintf("Experimental bin: %s\n\n", paste(exp_samples,  collapse = ", ")))

# --- Replicate assignment --------------------------------------------------
rep_guess_ctrl <- guess_replicates(ctrl_samples)
if (!is.null(rep_guess_ctrl)) {
  cat("Guessed replicate labels for control bin:\n")
  for (i in seq_along(ctrl_samples))
    cat(sprintf("  %s  ->  %s\n", ctrl_samples[i], rep_guess_ctrl[i]))
  if (toupper(trimws(readline("Correct? (Y/N): "))) != "Y")
    rep_guess_ctrl <- manual_assign_reps(ctrl_samples, "control")
} else {
  cat("Could not detect replicate numbers for control bin.\n")
  rep_guess_ctrl <- manual_assign_reps(ctrl_samples, "control")
}

rep_guess_exp <- guess_replicates(exp_samples)
if (!is.null(rep_guess_exp)) {
  cat("\nGuessed replicate labels for experimental bin:\n")
  for (i in seq_along(exp_samples))
    cat(sprintf("  %s  ->  %s\n", exp_samples[i], rep_guess_exp[i]))
  if (toupper(trimws(readline("Correct? (Y/N): "))) != "Y")
    rep_guess_exp <- manual_assign_reps(exp_samples, "experimental")
} else {
  cat("Could not detect replicate numbers for experimental bin.\n")
  rep_guess_exp <- manual_assign_reps(exp_samples, "experimental")
}

# --- Build col_data and DESeq2 object --------------------------------------
col_data <- data.frame(
  sample    = c(ctrl_samples, exp_samples),
  condition = c(rep("control", length(ctrl_samples)),
                rep("exp",     length(exp_samples))),
  replicate = c(rep_guess_ctrl, rep_guess_exp),
  row.names = c(ctrl_samples, exp_samples)
)
col_data$condition <- factor(col_data$condition, levels = c("control", "exp"))

cat("\nFinal sample assignments:\n")
print(col_data[, c("condition", "replicate")])
cat("\n")

if (toupper(trimws(readline("Proceed with these assignments? (Y/N): "))) != "Y")
  stop("Aborted by user. Re-run to reassign samples.")

counts_test <- counts_test[, rownames(col_data)]

dds <- DESeqDataSetFromMatrix(
  countData = counts_test,
  colData   = col_data,
  design    = ~ condition
)

write(c("", "Sample assignments:",
        capture.output(print(col_data[, c("condition", "replicate")]))),
      file = log_path, append = TRUE)

# =============================================================================
# 4. PRE-FILTERING
#
# Remove guides with very low counts before dispersion estimation.
# Guides with fewer than min_count counts in fewer than min_samples samples
# are excluded. These guides contribute little information and can
# destabilize dispersion estimation. The filter is applied to all guide types
# (test, pos, and neg).
# =============================================================================

min_count   <- 10  # adjust if needed: minimum counts per guide per sample
min_samples <- 3   # adjust if needed: minimum samples meeting min_count

keep           <- rowSums(counts(dds) >= min_count) >= min_samples
dropped_guides <- rownames(dds)[!keep]
dds            <- dds[keep, ]
cat("Test guides after low-count filter:", nrow(dds), "\n")

# The same filter is applied to the control (pos and neg) guides, which are
# analyzed separately in section 11.
ctrl_types   <- annotations$Guide_Type[annotations$Guide_Type %in% c("pos", "neg")]
ctrl_ids_all <- annotations$Guide_ID[annotations$Guide_Type %in% c("pos", "neg")]
ctrl_keep    <- rowSums(counts_raw[ctrl_ids_all, rownames(col_data), drop = FALSE] >=
                          min_count) >= min_samples
ctrl_dropped <- ctrl_ids_all[!ctrl_keep]
cat("Control guides removed by low-count filter:", length(ctrl_dropped), "\n\n")

write(c("",
        "Pre-filtering (low-count removal):",
        sprintf("  min_count:          %d", min_count),
        sprintf("  min_samples:        %d", min_samples),
        sprintf("  Test guides retained:    %d", nrow(dds)),
        sprintf("  Test guides dropped:     %d", length(dropped_guides)),
        if (length(dropped_guides) > 0)
          c("  Dropped test guide IDs:",
            paste0("    ", head(dropped_guides, 50)),
            if (length(dropped_guides) > 50)
              sprintf("    ... and %d more (full list omitted from log)",
                      length(dropped_guides) - 50)),
        sprintf("  Control guides dropped:  %d", length(ctrl_dropped)),
        if (length(ctrl_dropped) > 0)
          c("  Dropped control guide IDs:",
            paste0("    ", head(ctrl_dropped, 50)),
            if (length(ctrl_dropped) > 50)
              sprintf("    ... and %d more (full list omitted from log)",
                      length(ctrl_dropped) - 50))),
      file = log_path, append = TRUE)

# --- Guide accounting --------------------------------------------------------
# For each guide type: how many guides the annotations file lists, how many
# were absent from the counts file, how many were present but removed by the
# low-count filter, and how many were analyzed.
n_low_by_type <- c(test = length(dropped_guides),
                   pos  = sum(ctrl_types[!ctrl_keep] == "pos"),
                   neg  = sum(ctrl_types[!ctrl_keep] == "neg"))
acct <- data.frame(
  type      = guide_types,
  annotated = as.integer(n_annotated_by_type),
  missing   = as.integer(n_missing_by_type),
  low_count = as.integer(n_low_by_type[guide_types])
)
acct$analyzed <- acct$annotated - acct$missing - acct$low_count

acct_lines <- c(
  "",
  "Guide accounting (relative to the annotations file):",
  sprintf("  %-6s %10s %18s %20s %10s",
          "type", "annotated", "absent_from_counts", "removed_low_count", "analyzed"),
  sprintf("  %-6s %10d %18d %20d %10d",
          acct$type, acct$annotated, acct$missing, acct$low_count, acct$analyzed),
  sprintf("  Guides in the counts file without an annotation (excluded): %d",
          length(missing_ann)),
  if (length(missing_counts) > 0)
    c("  Annotated guides absent from the counts file:",
      paste0("    ", head(missing_counts, 50)),
      if (length(missing_counts) > 50)
        sprintf("    ... and %d more (full list omitted from log)",
                length(missing_counts) - 50))
)
write(acct_lines, file = log_path, append = TRUE)
cat(paste(acct_lines[-1], collapse = "\n"), "\n\n")

# =============================================================================
# 5. ITERATIVE SIGMA-CLIPPING NORMALIZATION
#
# Standard DESeq2 median-of-ratios normalization assumes most guides are
# unchanged between conditions. In noncoding tiling screens this assumption
# can fail: guides targeting functional elements may shift predominantly in
# one direction (depletion or enrichment, depending on screen modality and
# sort), biasing size factors estimated from all guides. Negative-control
# guides are an alternative reference, but they may not share the background
# behavior of targeting guides (for example, effects of DNA cutting).
#
# Size factors are therefore estimated from an empirically defined null set
# of test guides, identified by iterative sigma clipping:
#
#   Inner loop: Starting from all guides with valid LFC values, iteratively
#   exclude guides beyond ±k SD of the current mean and recompute until
#   membership converges (up to 50 inner iterations).
#
#   Outer loop: Use the converged null set as DESeq2 controlGenes to
#   reestimate size factors, rerun the Wald test, and reclip. Repeat until
#   null set membership is stable across passes (Jaccard > jaccard_thresh).
#
# The default clipping threshold is k = 2.0.
# =============================================================================

# --- Sigma clipping normalization function ---------------------------------

run_sigma_iterative <- function(dds,
                                k              = 2.0,
                                n_passes       = 5,
                                jaccard_thresh = 0.98,
                                min_null_n     = 100,
                                verbose        = TRUE) {

  # Dispersions are estimated once, under default size factors, and held
  # fixed across normalization passes.
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  null_idx_prev <- NULL
  pass_log      <- character(0)
  converged     <- FALSE      # set TRUE when the Jaccard threshold is reached
  last_jaccard  <- NA_real_   # most recent between-pass Jaccard overlap

  for (pass in seq_len(n_passes)) {

    dds <- nbinomWaldTest(dds)
    res <- results(dds, contrast = c("condition", "exp", "control"))

    lfc_vals <- res$log2FoldChange

    # Inner sigma clipping loop
    keep <- !is.na(lfc_vals)
    for (clip_iter in seq_len(50)) {
      mu       <- mean(lfc_vals[keep], na.rm = TRUE)
      sig      <- sd(lfc_vals[keep],   na.rm = TRUE)
      new_keep <- !is.na(lfc_vals) &
        lfc_vals >= (mu - k * sig) &
        lfc_vals <= (mu + k * sig)
      if (identical(new_keep, keep)) break
      keep <- new_keep
    }

    null_idx  <- which(keep)
    null_mean <- mean(lfc_vals[null_idx], na.rm = TRUE)
    null_sd   <- sd(lfc_vals[null_idx],   na.rm = TRUE)

    if (length(null_idx) < min_null_n) {
      warning(sprintf(
        "Pass %d: null set too small (%d guides). Consider increasing k.",
        pass, length(null_idx)
      ))
      pass_log <- c(pass_log, sprintf(
        "Pass %d: null set too small (%d) — stopped early.", pass, length(null_idx)
      ))
      break
    }

    pass_line <- sprintf(
      "Pass %d: null mean = %+.4f, null SD = %.4f, null n = %d (%.1f%%)",
      pass, null_mean, null_sd, length(null_idx),
      100 * length(null_idx) / nrow(dds)
    )

    if (verbose) cat(pass_line, "\n")
    pass_log <- c(pass_log, pass_line)

    # Convergence check
    if (!is.null(null_idx_prev)) {
      jaccard <- length(intersect(null_idx, null_idx_prev)) /
        length(union(null_idx, null_idx_prev))
      jaccard_line <- sprintf("  Jaccard overlap: %.3f", jaccard)
      if (verbose) cat(jaccard_line, "\n")
      pass_log <- c(pass_log, jaccard_line)
      last_jaccard <- jaccard
      if (jaccard >= jaccard_thresh) {
        converged <- TRUE
        if (verbose) cat("  Converged.\n")
        pass_log <- c(pass_log, "  Converged.")
        break
      }
    }

    dds <- estimateSizeFactors(dds, controlGenes = null_idx)
    null_idx_prev <- null_idx
  }

  # Warn if the loop ended without reaching the Jaccard threshold, either
  # because n_passes was exhausted or because the null set became too small.
  # In this case normalization has not fully stabilized; consider increasing
  # n_passes or k.
  if (!converged) {
    jaccard_txt <- if (is.na(last_jaccard)) "not computed" else sprintf("%.3f", last_jaccard)
    warning(sprintf(
      "Normalization did not converge (stopped after pass %d of %d; last Jaccard = %s). Consider increasing n_passes or k.",
      pass, n_passes, jaccard_txt
    ))
    pass_log <- c(pass_log, sprintf(
      "WARNING: did not converge (stopped after pass %d of %d; last Jaccard = %s).",
      pass, n_passes, jaccard_txt
    ))
  }

  # Run a final synchronized Wald test using the converged size factors so that
  # res and dds are always computed from the same normalization, regardless of
  # how the loop exited.
  dds <- nbinomWaldTest(dds)
  res <- results(dds, contrast = c("condition", "exp", "control"))

  list(
    dds      = dds,
    res      = res,
    null_idx = null_idx,
    null_mean = null_mean,
    null_sd   = null_sd,
    pass_log = pass_log
  )
}

# --- Normalization settings (adjust here if needed) -------------------------
k_clip         <- 2.0    # sigma clipping threshold (SD units)
n_passes       <- 5      # maximum outer normalization passes
jaccard_thresh <- 0.98   # null set overlap required for convergence
min_null_n     <- 100    # minimum null set size

k_txt <- format(k_clip, nsmall = 1)

# --- Run normalization -----------------------------------------------------
cat(sprintf("--- Iterative sigma clipping normalization (k = %s) ---\n", k_txt))
output_norm <- run_sigma_iterative(
  dds,
  k              = k_clip,
  n_passes       = n_passes,
  jaccard_thresh = jaccard_thresh,
  min_null_n     = min_null_n,
  verbose        = TRUE
)
cat("\n")

dds_final <- output_norm$dds
res_final  <- output_norm$res

# Append normalization settings and diagnostics to run log
write(c("", sprintf("Normalization strategy: sigma clipping (k = %s)", k_txt),
        sprintf("  k:                  %s", k_txt),
        sprintf("  n_passes:           %d", n_passes),
        sprintf("  jaccard_thresh:     %s", format(jaccard_thresh)),
        sprintf("  min_null_n:         %d", min_null_n),
        "--- Iterative sigma clipping normalization ---",
        output_norm$pass_log),
      file = log_path, append = TRUE)

# =============================================================================
# 6. SHRUNKEN LFC (apeglm)
#
# lfcShrink applies empirical Bayes shrinkage to MLE LFC estimates using the
# converged sigma-clipping size factors, placing shrunken values on the same
# scale as the Wald test results. Shrunken LFC is preferred for ranking and
# visualization; significance calls are derived from the Wald statistic and
# are unaffected by shrinkage.
# =============================================================================

cat("--- LFC shrinkage (apeglm) ---\n")
apeglm_coef <- resultsNames(dds_final)[2]
cat("Coefficient:", apeglm_coef, "\n")

res_shrunk <- lfcShrink(
  dds_final,
  coef = apeglm_coef,
  type = "apeglm"
)
cat("\n")

# =============================================================================
# 7. NORMALIZED COUNTS AND PER-REPLICATE LFC HELPER
#
# For any DESeq2 object, computes:
#   norm_*    — normalized count per sample (raw count / size factor)
#   L2FC_Rep* — per-replicate LFC: log2((norm_exp+1) / (norm_ctrl+1))
#               pseudocount of +1 handles zero counts; replicates matched
#               by label
# =============================================================================

get_counts_and_replfc <- function(dds) {

  norm_counts <- counts(dds, normalized = TRUE)
  cd          <- as.data.frame(colData(dds))

  ctrl_cols        <- rownames(cd)[cd$condition == "control"]
  exp_cols         <- rownames(cd)[cd$condition == "exp"]
  matched          <- match(cd[ctrl_cols, "replicate"], cd[exp_cols, "replicate"])
  exp_cols_matched <- exp_cols[matched]

  if (any(is.na(matched)))
    warning("Could not match all replicates by label — check replicate assignments")

  norm_df <- norm_counts %>%
    as.data.frame() %>%
    rownames_to_column("guide")
  colnames(norm_df)[-1] <- paste0("norm_", colnames(norm_df)[-1])

  rep_lfc_df <- map2_dfc(ctrl_cols, exp_cols_matched, function(ctrl, ex) {
    rep_label <- cd[ctrl, "replicate"]
    tibble(!!paste0("L2FC_", rep_label) :=
             log2((norm_counts[, ex] + 1) / (norm_counts[, ctrl] + 1)))
  }) %>%
    mutate(guide = rownames(norm_counts)) %>%
    select(guide, everything())

  left_join(norm_df, rep_lfc_df, by = "guide")
}

# =============================================================================
# 8. COMPILE TEST GUIDE RESULTS TABLE
# =============================================================================

# --- One-sided p-values from Wald statistic --------------------------------
# The Wald statistic is asymptotically standard normal under the null,
# allowing derivation of one-sided p-values without additional assumptions.
# BH controls FDR under independence or positive dependence.
#
# Note on outliers: DESeq2's Cook's distance outlier check (and its
# independent filtering) set only DESeq2's own two-sided pvalue/padj to NA;
# the Wald statistic is reported for every guide. Because the one-sided
# p-values here are computed from the Wald statistic, Cook's distance outlier
# filtering is not applied: every test guide that passes the low-count filter
# is tested and included in the BH correction.

stat_vec <- res_final$stat

p_deplete         <- pnorm(stat_vec, lower.tail = TRUE)
padj_depletion_BH <- p.adjust(p_deplete, method = "BH")

p_enrich           <- pnorm(stat_vec, lower.tail = FALSE)
padj_enrichment_BH <- p.adjust(p_enrich, method = "BH")

# --- Normalized counts and per-replicate LFC -------------------------------
counts_replfc_test <- get_counts_and_replfc(dds_final)

# --- Assemble results table ------------------------------------------------
results_final <- res_final %>%
  as.data.frame() %>%
  rownames_to_column("guide") %>%
  left_join(
    res_shrunk %>%
      as.data.frame() %>%
      rownames_to_column("guide") %>%
      select(guide, log2FC_shrink = log2FoldChange, lfcSE_shrunk = lfcSE),
    by = "guide"
  ) %>%
  mutate(
    p_deplete          = p_deplete,
    padj_depletion_BH  = padj_depletion_BH,
    p_enrich           = p_enrich,
    padj_enrichment_BH = padj_enrichment_BH
  ) %>%
  left_join(counts_replfc_test, by = "guide") %>%
  left_join(annotations %>% rename(guide = Guide_ID), by = "guide") %>%
  mutate(
    coord        = Coordinate,
    type         = Guide_Type,
    in_sigma_null = guide %in% rownames(dds_final)[output_norm$null_idx]
  ) %>%
  select(guide, coord,
         type,
         any_of(setdiff(colnames(annotations), c("Guide_ID", "Guide_Type", "Coordinate"))),
         baseMean,
         log2FoldChange, lfcSE, stat, pvalue, padj,
         log2FC_shrink, lfcSE_shrunk,
         p_deplete,  padj_depletion_BH,
         p_enrich,   padj_enrichment_BH,
         in_sigma_null,
         starts_with("norm_"),
         starts_with("L2FC_")) %>%
  arrange(log2FoldChange)

# =============================================================================
# 9. RESULTS SUMMARY
# =============================================================================

count_calls <- function(df, padj_col, threshold) {
  sum(df[[padj_col]] < threshold, na.rm = TRUE)
}

summary_lines <- c(
  "=============================================================================",
  "RESULTS SUMMARY",
  "=============================================================================",
  sprintf("Normalization strategy: sigma clipping (k = %s)", k_txt),
  sprintf("Total test guides tested: %d", nrow(results_final)),
  sprintf("Null set size:            %d (%.1f%%)",
          sum(results_final$in_sigma_null, na.rm = TRUE),
          100 * mean(results_final$in_sigma_null, na.rm = TRUE)),
  "",
  "--- One-sided depletion test ---",
  sprintf("  BH < 0.05:  %d", count_calls(results_final, "padj_depletion_BH", 0.05)),
  sprintf("  BH < 0.01:  %d", count_calls(results_final, "padj_depletion_BH", 0.01)),
  "",
  "--- One-sided enrichment test ---",
  sprintf("  BH < 0.05:  %d", count_calls(results_final, "padj_enrichment_BH", 0.05)),
  sprintf("  BH < 0.01:  %d", count_calls(results_final, "padj_enrichment_BH", 0.01))
)

cat(paste(summary_lines, collapse = "\n"), "\n\n")
write(c("", paste(summary_lines, collapse = "\n")),
      file = log_path, append = TRUE)

# =============================================================================
# 10. PCA
# =============================================================================

# --- PCA on VST-normalized counts ------------------------------------------
# Runs on all post-filter test guides. The VST stabilizes variance across the
# count range so that high-count guides do not dominate the PCA axes.
# Samples should cluster by condition (control vs exp) along PC1 if the
# biological signal is consistent across replicates. Poor clustering or
# swapped replicates would be apparent in the PC scores written to QC_PCA.tsv.

vst_data <- vst(dds_final, blind = FALSE)

# Write PCA coordinates to a single TSV for plotting in external software.
# Uses the 500 highest-variance guides (the DESeq2 plotPCA default).
# Columns: sample, condition, replicate, then one column per component,
# with the percent variance explained in the header, e.g. "PC1 (89.24%)".

pca_input  <- assay(vst_data)
top_var    <- head(order(matrixStats::rowVars(pca_input), decreasing = TRUE), 500)
pca_result <- prcomp(t(pca_input[top_var, ]))

pct_var <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 2)

pca_scores <- pca_result$x
colnames(pca_scores) <- sprintf("%s (%.2f%%)", colnames(pca_scores), pct_var)

pca_out <- col_data[rownames(pca_scores), c("condition", "replicate")] %>%
  rownames_to_column("sample") %>%
  bind_cols(as.data.frame(pca_scores, check.names = FALSE))

pca_path <- file.path(output_dir, "QC_PCA.tsv")
write.table(pca_out, file = pca_path,
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("PCA coordinates written to: %s\n", pca_path))
cat(sprintf("  PC1: %.1f%% variance  PC2: %.1f%% variance\n",
            pct_var[1], pct_var[2]))

# =============================================================================
# 11. CONTROL GUIDE RESULTS
#
# Positive (pos) and negative (neg) control guides are processed using the
# size factors from the converged test guide normalization, placing their
# LFC values on the same scale as test guide results. Control guides removed
# by the low-count filter (section 4) are excluded.
#
#   pos: positive controls — analyzed for validation and quality control
#   neg: negative controls — analyzed for validation and quality control
# =============================================================================

process_control_guides <- function(guide_type, counts_raw, annotations,
                                   col_data, size_factors, exclude_ids) {

  cat(sprintf("\n--- Processing %s guides ---\n", guide_type))

  ctrl_guides <- annotations %>%
    filter(Guide_Type == guide_type) %>%
    pull(Guide_ID)
  n_found     <- length(ctrl_guides)
  ctrl_guides <- setdiff(ctrl_guides, exclude_ids)   # low-count filter (section 4)

  counts_ctrl <- counts_raw[ctrl_guides, rownames(col_data), drop = FALSE]
  cat(sprintf("  %s guides found: %d (%d removed by low-count filter, %d analyzed)\n",
              guide_type, n_found, n_found - nrow(counts_ctrl), nrow(counts_ctrl)))

  dds_ctrl <- DESeqDataSetFromMatrix(
    countData = counts_ctrl, colData = col_data, design = ~ condition
  )
  sizeFactors(dds_ctrl) <- size_factors
  dds_ctrl <- estimateDispersions(dds_ctrl)
  dds_ctrl <- nbinomWaldTest(dds_ctrl)
  res_ctrl  <- results(dds_ctrl, contrast = c("condition", "exp", "control"))

  res_ctrl_shrunk <- lfcShrink(
    dds_ctrl, coef = resultsNames(dds_ctrl)[2], type = "apeglm"
  )

  stat_vec           <- res_ctrl$stat
  p_deplete          <- pnorm(stat_vec, lower.tail = TRUE)
  padj_depletion_BH  <- p.adjust(p_deplete,  method = "BH")
  p_enrich           <- pnorm(stat_vec, lower.tail = FALSE)
  padj_enrichment_BH <- p.adjust(p_enrich, method = "BH")
  # Note: FDR correction here uses only the control guide set as the multiple
  # testing universe, not the full test guide set. Control guide adjusted
  # p-values are therefore not directly comparable to test guide padj values
  # and should be used for QC validation only.

  counts_replfc_ctrl <- get_counts_and_replfc(dds_ctrl)

  results_ctrl <- res_ctrl %>%
    as.data.frame() %>%
    rownames_to_column("guide") %>%
    left_join(
      res_ctrl_shrunk %>%
        as.data.frame() %>%
        rownames_to_column("guide") %>%
        select(guide, log2FC_shrink = log2FoldChange, lfcSE_shrunk = lfcSE),
      by = "guide"
    ) %>%
    mutate(
      p_deplete          = p_deplete,
      padj_depletion_BH  = padj_depletion_BH,
      p_enrich           = p_enrich,
      padj_enrichment_BH = padj_enrichment_BH
    ) %>%
    left_join(counts_replfc_ctrl, by = "guide") %>%
    left_join(annotations %>% rename(guide = Guide_ID), by = "guide") %>%
    mutate(coord = Coordinate, type = Guide_Type) %>%
    select(guide, coord, type,
           any_of(setdiff(colnames(annotations), c("Guide_ID", "Guide_Type", "Coordinate"))),
           baseMean, log2FoldChange, lfcSE, stat, pvalue, padj,
           log2FC_shrink, lfcSE_shrunk,
           p_deplete,  padj_depletion_BH,
           p_enrich,   padj_enrichment_BH,
           starts_with("norm_"), starts_with("L2FC_")) %>%
    arrange(log2FoldChange)

  cat(sprintf("  Median LFC:                            %+.3f\n",
              median(results_ctrl$log2FoldChange, na.rm = TRUE)))
  cat(sprintf("  Depleted (padj_depletion_BH < 0.05):   %d\n",
              sum(results_ctrl$padj_depletion_BH < 0.05, na.rm = TRUE)))
  cat(sprintf("  Enriched (padj_enrichment_BH < 0.05):  %d\n",
              sum(results_ctrl$padj_enrichment_BH < 0.05, na.rm = TRUE)))

  invisible(results_ctrl)
}

results_pos <- process_control_guides("pos", counts_raw, annotations,
                                      col_data, sizeFactors(dds_final), ctrl_dropped)
results_neg <- process_control_guides("neg", counts_raw, annotations,
                                      col_data, sizeFactors(dds_final), ctrl_dropped)

# =============================================================================
# 12. WRITE RESULTS TABLES
#
# Results are written to two tables with the same column layout:
#   DESeq2_results_test.tsv      — test guides (input to the window analysis)
#   DESeq2_results_controls.tsv  — positive and negative control guides
# Rows follow the order of the counts file. Adjusted p-values are
# Benjamini-Hochberg. Control guide adjusted p-values are corrected within
# their own guide type (see section 11) and are intended for QC only, so they
# are not comparable to test guide adjusted p-values. 'DESeq2 Null Group'
# (test table only) is TRUE if the guide was in the converged null set.
# =============================================================================

null_ids <- rownames(dds_final)[output_norm$null_idx]

# Match the annotation column settings to the file's column names without
# regard to capitalization; the file's own spelling is used in the output.
match_ann_cols <- function(wanted) {
  found <- colnames(annotations)[match(tolower(wanted), tolower(colnames(annotations)))]
  found[!is.na(found)]
}
missing_ann_cols <- c(annotation_cols_front, annotation_cols_back)[
  !tolower(c(annotation_cols_front, annotation_cols_back)) %in% tolower(colnames(annotations))]
if (length(missing_ann_cols) > 0)
  warning(sprintf("Annotation column(s) not found and omitted from output: %s",
                  paste(missing_ann_cols, collapse = ", ")))
annotation_cols_front <- match_ann_cols(annotation_cols_front)
annotation_cols_back  <- match_ann_cols(annotation_cols_back)

results_table <- bind_rows(results_final, results_pos, results_neg) %>%
  mutate(`DESeq2 Null Group` = ifelse(type == "test", guide %in% null_ids, NA)) %>%
  select(GuideID          = guide,
         type,
         Coordinate       = coord,
         any_of(annotation_cols_front),
         `Raw L2FC`       = log2FoldChange,
         `Raw L2FC SE`    = lfcSE,
         stat,
         `apeglm L2FC`    = log2FC_shrink,
         `apeglm L2FC SE` = lfcSE_shrunk,
         `p deplete`      = p_deplete,
         `adj-p deplete`  = padj_depletion_BH,
         `p enrich`       = p_enrich,
         `adj-p enrich`   = padj_enrichment_BH,
         baseMean,
         starts_with("norm_"),
         starts_with("L2FC_"),
         any_of(annotation_cols_back),
         `DESeq2 Null Group`) %>%
  rename_with(~ paste0("norm-cnt ", sub("^norm_", "", .x)),
              starts_with("norm_")) %>%
  rename_with(~ sub("^L2FC_", "L2FC ", .x), starts_with("L2FC_")) %>%
  arrange(match(GuideID, rownames(counts_raw)))

results_test     <- results_table %>% filter(type == "test")
results_controls <- results_table %>% filter(type != "test") %>%
  select(-`DESeq2 Null Group`)

test_path     <- file.path(output_dir, "DESeq2_results_test.tsv")
controls_path <- file.path(output_dir, "DESeq2_results_controls.tsv")
write.table(results_test, file = test_path,
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(results_controls, file = controls_path,
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")

cat(sprintf("\nTest guide results written to:    %s (%d guides)\n",
            test_path, nrow(results_test)))
cat(sprintf("Control guide results written to: %s (%d guides; pos: %d, neg: %d)\n",
            controls_path, nrow(results_controls),
            sum(results_controls$type == "pos"),
            sum(results_controls$type == "neg")))

cat("\nDone.\n")
