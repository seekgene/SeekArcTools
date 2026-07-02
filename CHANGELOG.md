# Changelog
## 2026-06-25

- Modifies the way cell identification results from running cell_identify.R are processed in the astep3.py script.

## 2026-06-22

### Removed

- Removed deprecated `--filter_gex_umi` option and related GEX UMI post-filtering logic from SeekARC ARC analysis.

## 2026-06-12

### Added

- Added optional `pairendAlignment` support for SeekARC RNA-only analysis.

  - Default RNA-only alignment remains single-end R2 alignment.
  - When enabled, step1 output R1/R2 FASTQ files are used for STAR paired-end alignment.

- Added synchronized R1/R2 filtering before paired-end alignment.

  - Read pairs with R1 shorter than `pairendMinR1Len` are removed together to avoid STAR failure caused by empty or too-short R1 after barcode/UMI trimming.
  - `qualimap`, `featureCounts`, and RNA UMI counting are updated to handle paired-end mode.

### Updated

* Updated `seekarctools rna run` to execute the complete RNA-only workflow.

  * After RNA counting and cell calling, `run` now automatically generates Seurat RDS and RNA-only GEX summary CSV.
  * Added `memory` parameter support in `rna run` for Seurat RDS generation.

## 2026-06-11

### Added

- Added SeekARC RNA-only software support.

  - Added `DD_AG_RNA` chemistry configuration for SeekARC RNA-only data.
  - Added the `seekarctools rna` command group for GEX-only analysis.
  - Added RNA-only step-wise commands: `step1`, `step2`, `step3`, `rds`, `report`, and `run`.
  - Added STAR loose-alignment parameter passthrough with `scoremin` and `matchnmin`.
  - Added RNA-only cell calling from `raw_feature_bc_matrix` to `filtered_feature_bc_matrix`.
  - Added RNA-only Seurat RDS generation.
  - Added RNA-only GEX summary CSV generation.

- Added a compatibility helper for ARC GEX cell calling.

  - Added `cell_calling()` in `seekarctools.arc.estep3`.
  - This helper calls `utils/cell_identify.R` to generate `filtered_feature_bc_matrix` from `raw_feature_bc_matrix`.
  - The existing `arc estep3` counting behavior remains unchanged.

### Fixed

- Fixed the missing `cell_calling()` definition in `seekarctools.arc.estep3`, which could cause errors when downstream code imports this helper.


## 2026-06-03

- Enhanced Cell Calling: Introduced `strict_cell_calling` parameter, integrating the EmptyDrops algorithm for improved accuracy.
- Flexible Cell Limits: Replaced the hardcoded cell count limit with a customizable `max_cells` parameter.
- Post-processing: Added `filter_gex_umi` for downstream GEX UMI filtering.
- Reproducibility: Fixed the random seed for K-Means clustering to ensure consistent results across runs.
- Parameter Renaming: Updated `mito_filter` to `keep_mito` for better clarity.
- CLI Optimization: Streamlined the `seekarctools utils help` menu by hiding low-level reference genome commands (e.g., `mkref`, `gtffilter`, `splitref`).

## 2026-05-26

- Data Integration: Mitochondrial metrics are now exported to `summary.csv` and synchronized with cloud platforms and databases.
- GTF Parsing Engine Upgrade (`count_link`):

  - Refactored Logic: Re-engineered the underlying GTF reading mechanism specifically for `astep4` linkages.
  - Performance & Robustness: Optimized parsing significantly reduces peak memory consumption and enhances compatibility with non-standard or complex GTF formats (e.g., those with redundant tags), preventing parsing failures or feature loss.
