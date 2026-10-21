# SPANA: Saturating Perturbation Analysis of Noncoding Activity

Analysis code accompanying:

> Vyas HS, Leix K, Burugula BB, Paye SM, Raut C, Speliotes EK, Emmer BT. Multimodal functional characterization of the LDLR noncoding landscape. Manuscript in preparation.

## Overview

SPANA consists of three R scripts that analyze pooled noncoding CRISPR tiling screens with two FACS-sorted populations. SPANA-null runs first; SPANA-scan and SPANA-cluster each use its test-guide results. They were used for the CRISPRn, CRISPRi, and CRISPRa screens in the paper.

1. **SPANA-null** (`SPANA-null.R`) identifies an empirical null set of test guides, uses it to normalize the screen, and then tests each guide for enrichment or depletion with DESeq2.
2. **SPANA-scan** (`SPANA-scan.R`) averages guide-level effects across sliding genomic windows and calls peaks.
3. **SPANA-cluster** (`SPANA-cluster.R`) tests whether significantly depleted or enriched guides cluster along the locus more than expected by chance.

## Requirements

- R (the paper used R 4.4.3)
- CRAN: `tidyverse`
- Bioconductor: `DESeq2`, `apeglm` (`matrixStats` is installed automatically with DESeq2)

```r
install.packages(c("tidyverse", "BiocManager"))
BiocManager::install(c("DESeq2", "apeglm"))
```

SPANA-scan and SPANA-cluster use only base R. Package versions are recorded in each run log.

## Step 1: SPANA-null (null-set normalization and guide-level analysis)

**Normalization.** Standard median-of-ratios normalization assumes most guides are unchanged between conditions. In noncoding tiling screens, guides targeting functional sequences may shift predominantly in one direction, biasing size factors estimated from all guides. Size factors are therefore estimated from an empirically defined null set of test guides, identified by iterative sigma clipping of the log2 fold change (LFC) distribution. Guides beyond ±*k* SD of the mean (default *k* = 2) are excluded, and the mean and SD are recomputed until the set is stable. Size factors are then re-estimated from that set, and the process repeats until null-set membership is stable between passes (Jaccard overlap ≥ 0.98).

**Significance testing.** One-sided p-values for depletion and enrichment are derived from the DESeq2 Wald statistic and adjusted using the Benjamini–Hochberg (BH) method. Cook's distance outlier filtering is not applied to the one-sided tests. Shrunken LFC estimates are computed with apeglm.

**Guide types.** Guides labeled `test` are used for normalization and hit calling. Positive (`pos`) and negative (`neg`) controls are not used for normalization. They are analyzed separately, for validation and quality control, using the size factors from the test guides. Their adjusted p-values are corrected within each control type.

**Filtering.** Before analysis, guides are removed unless at least 3 samples have 10 or more counts. This applies to all guide types.

### Input files

Both files are tab-delimited and must be in the working directory.

**1. Counts file.** The first column holds guide IDs, and every other column is read as a sample, so the file should contain no other columns.

| GuideID | CRISPRi LDL-low Rep1 | ... | CRISPRi LDL-high Rep4 |
|---|---|---|---|
| LDLR_11039456+ | 973 | ... | 349 |

**2. Annotations file.** One row per guide. Column names are matched without regard to capitalization.

| Column | Required | Description |
|---|---|---|
| `guide_ID` | Yes | Must match the guide IDs in the counts file |
| `guide_type` | Yes | `test`, `pos`, or `neg` (capitalization and surrounding spaces are ignored; any other value stops the script) |
| `coordinate` | Yes | Genomic cut-site coordinate; may be blank or `NA` for guides without one |
| Any others | No | Can be carried through to the output; see "Settings" below |

Numbers written with thousands separators (for example `11,039,456`) are read as plain numbers.

Guides in the counts file without an annotation, and annotated guides missing from the counts file, are excluded with a warning. Both are reported in the run log.

### Running

```r
setwd("path/to/directory/with/input/files")
source("SPANA-null.R")
```

The script asks you to:

1. Choose the counts file and the annotations file from a numbered list
2. Confirm the two sample groups. These are detected from sample names (for example `high`/`low`), or you can enter them manually.
3. Choose which group is the control (reference). Positive LFC means enriched in the other group.
4. Confirm the replicate labels for each group, which are used to pair samples for the per-replicate LFC columns
5. Confirm the final sample assignments

### Settings

**Output settings** (top of the script). `annotation_cols_front` and `annotation_cols_back` name the annotation columns to include in the results tables, placed after `Coordinate` and after the `L2FC <rep>` columns, respectively. Names are matched without regard to capitalization. The defaults are `group` and `LDLR CDS proximity`.

**Analysis settings.**

| Setting | Default | Location | Description |
|---|---|---|---|
| `k_clip` | 2.0 | Section 5 | Sigma-clipping threshold (SD units) |
| `n_passes` | 5 | Section 5 | Maximum normalization passes |
| `jaccard_thresh` | 0.98 | Section 5 | Null-set overlap required for convergence |
| `min_null_n` | 100 | Section 5 | Minimum null-set size |
| `min_count` | 10 | Section 4 | Minimum count per guide per sample |
| `min_samples` | 3 | Section 4 | Minimum number of samples meeting `min_count` |

### Outputs

Outputs are written to a timestamped folder (`YYYY_MM_DD_HHMM_DESeq2_output/`).

| File | Contents |
|---|---|
| `DESeq2_results_test.tsv` | Results for test guides; the input for Steps 2 and 3 |
| `DESeq2_results_controls.tsv` | Results for positive and negative control guides |
| `QC_PCA.tsv` | PCA of VST-normalized test-guide counts: PC scores per sample, with percent variance explained in each column header |
| `run_log.txt` | Script version, input file checksums, package versions, sample assignments, guide accounting, normalization settings and passes, and a results summary |

**Results table columns** (both results files)

| Column | Description |
|---|---|
| `GuideID`, `type`, `Coordinate` | Guide ID, guide type, and cut-site coordinate |
| `Raw L2FC`, `Raw L2FC SE` | Maximum-likelihood LFC and its standard error |
| `stat` | Wald statistic |
| `apeglm L2FC`, `apeglm L2FC SE` | apeglm-shrunken LFC and its posterior SD |
| `p deplete`, `adj-p deplete` | One-sided depletion p-value and its BH-adjusted value |
| `p enrich`, `adj-p enrich` | One-sided enrichment p-value and its BH-adjusted value |
| `baseMean` | Mean normalized count across samples |
| `norm-cnt <sample>` | Normalized count per sample, named as in the counts file |
| `L2FC <rep>` | Per-replicate LFC: log2((normalized experimental + 1) / (normalized control + 1)) |
| `DESeq2 Null Group` | Test file only: TRUE if the guide was in the converged null set |

Results may differ in the last decimal places between DESeq2/apeglm versions or with a different guide order in the input files.

## Step 2: SPANA-scan (sliding-window analysis)

For each test guide, the script averages the apeglm-shrunken LFC of all test guides with cut sites within ±49 bp (a 100-bp window). Each window mean is compared with a null distribution built by randomly sampling the same number of test-guide LFCs, giving empirical one-sided p-values for depletion and enrichment, which are BH-adjusted. Overlapping significant windows are merged into peaks. Control guides are not used.

### Running

Set the working directory to the Step 1 output folder, which contains `DESeq2_results_test.tsv`:

```r
setwd("path/to/YYYY_MM_DD_HHMM_DESeq2_output")
source("SPANA-scan.R")
```

The script asks you to confirm the input file and enter an output prefix, the number of simulations, the peak-calling FDR threshold, the chromosome of the screened locus, and the random seed. Defaults are offered for all except the prefix and chromosome. To run without prompts (for example on a computing cluster), fill in the parameters at the top of the script and set `INTERACTIVE_MODE <- FALSE`.

The simulations run in parallel on Linux and macOS (using one fewer core than are available or allocated) and on a single core on Windows. Results for a given seed are the same regardless of the number of cores.

### Settings

| Setting | Default | Description |
|---|---|---|
| `WINDOW_BP` | 100 | Window width (bp) |
| `N_SIM` | 100,000 | Simulations per window size |
| `FDR_THRS` | 0.01 | BH-adjusted p-value threshold for peak calling |
| `RANDOM_SEED` | 42 | Random seed |

### Outputs

Outputs are written to a timestamped subfolder (`YYYY_MMDD_HHMM_window_output/`).

| File | Contents |
|---|---|
| `<prefix>_window_results_LFCshr.tsv` | One row per window (centered on each test guide): `GuideID`, `Coordinate`, number of guides in the window, window mean apeglm LFC, and depletion and enrichment p-values with BH-adjusted values |
| `<prefix>_LFCshr_{depletion,enrichment}_BH_FDR001_peaks.bed` | Peaks at the FDR threshold (file names reflect the threshold used). Each window is written as a 100-bp interval centered on its guide, 1 bp wider than the span of cut sites it can include. |
| `<prefix>_window_mean_LFCshr.wig` and `<prefix>_window_empirical_p{Deplete,Enrich}_LFCshr.wig` | Genome browser tracks of window mean LFC and −log10 p-values |
| `<prefix>_run_parameters.txt` | Parameters, guide counts, and R session information |

## Step 3: SPANA-cluster (nearest-neighbor clustering test)

Test guides are called hits if their BH-adjusted p-value is below the FDR threshold. The test statistic is the median distance from each hit to its nearest neighboring hit. The null distribution is built by repeatedly drawing the same number of guides at random from all test guides, keeping their actual positions. The one-sided empirical p-value is (number of permutations with a median distance ≤ observed + 1) / (number of permutations + 1); when no permutation reaches the observed value, this is the smallest p-value the permutations can resolve. Depleted and enriched hits are tested separately, and a direction with fewer than 2 hits is skipped. Guides are sorted by coordinate before permuting, so results for a given seed do not depend on the row order of the input.

### Running

Set the working directory to the Step 1 output folder, which contains `DESeq2_results_test.tsv`:

```r
setwd("path/to/YYYY_MM_DD_HHMM_DESeq2_output")
source("SPANA-cluster.R")
```

The script asks you to confirm the input file and enter an output prefix, the FDR threshold, the number of permutations, and the random seed. Defaults are offered for all except the prefix. To run without prompts, fill in the parameters at the top of the script and set `INTERACTIVE_MODE <- FALSE`.

### Settings

| Setting | Default | Description |
|---|---|---|
| `FDR_THR` | 0.01 | Hits: BH-adjusted p-value below this threshold |
| `N_PERM` | 100,000 | Permutations per direction |
| `RANDOM_SEED` | 1 | Random seed |

### Outputs

Outputs are written to a timestamped subfolder (`YYYY_MMDD_HHMM_cluster_output/`).

| File | Contents |
|---|---|
| `<prefix>_nn_null_distributions.tsv` | Permutation null values of the median nearest-neighbor distance (bp), one column per direction |
| `<prefix>_nn_summary.tsv` | Per direction: number of hits and guides tested, observed median nearest-neighbor distance (bp), p-value, and number of permutations |
| `<prefix>_run_parameters.txt` | Parameters, guide counts, results, and R session information |

## Reproducing the paper's results

- Input data: guide annotations are in Supplementary Table 1, and read counts are in Supplementary Tables 2 (CRISPRn), 3 (CRISPRi), and 4 (CRISPRa), provided with the paper and the bioRxiv preprint. Save the annotation table and each screen's read-count sheet as tab-delimited text files before running.
- Step 1: choose the LDL-low samples as the control group in each screen.
- Step 2: chromosome 19, seed 42, FDR 0.01. The paper used 1,000,000,000 simulations per window size, which takes considerably longer than the default. Because the null is simulated, windows with adjusted p-values close to the threshold may be called differently with other settings.
- Step 3: the defaults (FDR 0.01, 100,000 permutations, seed 1).

## Citation

If you use this code, please cite the paper above and the archived software (see `CITATION.cff`, or use "Cite this repository" on GitHub).

## Funding

This research was supported by the National Institutes of Health R01-HL167733 (BTE), R01-HL171013 (BTE), the A. Alfred Taubman Medical Research Institute (BTE), the American Heart Association 24POST1198445 (HSV), and the Michigan Pioneer Fellows Program (HSV). The funders had no role in study design, data collection and analysis, decision to publish, or preparation of the manuscript.

## License

Copyright (C) 2026 The Regents of the University of Michigan
Emmer Lab

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See `LICENSE` for details, and `COPYRIGHT.txt` for the licenses of the third-party R packages SPANA uses.

## Contact

Brian T. Emmer (bemmer@med.umich.edu), Emmer Lab, University of Michigan
