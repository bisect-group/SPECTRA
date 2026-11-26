# BuildCancerModels — CCLE Context-Specific Cancer Metabolic Models & Flux Analysis

## Requirements
1. Cancer Cell Line Encyclopedia (CCLE) transcriptomic expression profiles for 1,479 cell lines (`DepMap Public 23Q4`, https://depmap.org/portal/).
2. MATLAB with the COBRA Toolbox (`fastCoreWeighted`), `SPECTRA` (`spectraME`), `SWIFTCORE`, and Gurobi (`feasTol = 1e-8`).
3. Python 3.10+ (`numpy`, `scipy`, `pandas`, `matplotlib`, `seaborn`) for statistical analyses, RHMC MCMC flux-sampling corroboration, and Figure 5 plotting.

---

## Directory Structure

### 1. `DepMapModels_Recon/` — Reconstructing 1,479 CCLE Context-Specific Models
- `PrepareDepMapData.m`: Extracts CCLE transcriptomic expression profiles into `DepMapData.mat`.
- `modGenesIDConv.m`: Maps Recon3D+ gene IDs (`UpdatedRecon3D.mat`) to gene symbols (`consRecon3DGeneSymbol.mat`).
- `GiniReactionImportance.m`: Computes `LocalGini` core reaction sets (`RxnImp >= 1`) and non-core penalty weights (`1 - RxnImp`).
- `BuildDepMapModels.m`: Reconstructs context-specific models across all 1,479 CCLE cell lines using weighted Fastcore (`fastCoreWeighted`), `swiftcore`, and `spectraME` (`minNetLP`).
- `run_fastCoreWeighted_array.m`: Parallel chunked SLURM runner for `fastCoreWeighted` across the 1,479 CCLE cell lines.
- `run_sp_new_array.m` & `run_timing_1479_array.m`: Parallel chunked SLURM runners for `spectraME` (`minNetLP`) model reconstruction and back-to-back 3-method extraction runtime benchmarking across the 1,479 CCLE cell lines.
- `run_fva20_array.m`: Parallel SLURM runner for unconstrained Flux Variability Analysis (`optPercentage = 0`) across the 20 paired cancer and non-cancerous `SPECTRA_ME` models.

### 2. `coverage_core/` — 135 Cancer Hallmark Recovery, Model Size `|N|`, & Specificity Metrics
- `cancerrxnscoverage.m` & `core_coverage_find.m`: Computes the proportion of 135 curated cancer hallmark reactions recovered in each extracted CSM (`coverage_spectra.mat`, `coverage_fastcore.mat` / `coverage_fastcore_weighted.mat`, `coverage_swiftcore.mat`).
- `ccle_1479_size_and_hallmark_metrics.csv`: Per-cell-line core size `|C|`, extracted model size `|N|`, non-core additions `|N \ C|`, hallmark recall, non-core hypergeometric fold-enrichment & Fisher odds ratios, negative-control pathway fold-enrichment, and back-to-back per-model runtimes across all 1,479 CCLE lines.
- `ccle_1479_wilcoxon_stats.csv`: Unpaired two-sided Mann-Whitney U and paired Wilcoxon signed-rank statistics across all 1,479 CCLE CSMs.
- `analyze_ccle_1479_and_fig5c.py`: Generates the 3-panel Figure 5c (`Visualisations/cancer_pareto_speed_coverage.pdf`: hallmark recovery, extracted model size, and extraction runtime).

### 3. `fva/` — Flux Variability Analysis Bounds (20 Paired Cancer vs. Normal CSMs)
- Contains `minFlux` and `maxFlux` FVA bounds (`ACH-*.xlsx`) for the 10 cancer and 10 tissue-matched non-cancerous `SPECTRA_ME` models listed in `flux_results/models_name.m`.

### 4. `flux_results/` — Flux Span Ratio (FSR), Hypergeometric FEA, & Figure 5d
- `models_name.m`: Lists the 10 cancer and 10 non-cancerous CCLE `ACH-...` identifiers.
- `readfva.m`: Aggregates FVA Excel bounds into `cells_recon_2025.mat`.
- `fsranalysis_ccle_sprint.m`: Identifies upregulated (`maxFlux_cancer > maxFlux_normal` & `FSR >= 2`) and downregulated (`maxFlux_cancer < maxFlux_normal` & `FSR <= 0.5`) reactions and runs one-sided upper-tail hypergeometric Flux Enrichment Analysis (`BH FDR < 0.05`).
- `fsr_up_2025_checked_2.xlsx`, `fsr_down_2025_checked_2.xlsx`, `fea_2025_checked_2.xlsx`: Multi-sheet workbooks corresponding to Supplementary Data S4, S5, and S6 (including FVA bounds and RHMC MCMC corroboration statistics).
- `step4_rxn_fsr_mcmc_concordance.csv` & `step4_new_fsr_pathways_with_mcmc_corroboration.csv`: Reaction- and pathway-level concordance between FSR/FEA and RHMC MCMC flux sampling.
- `run_fsr_fea_mcmc_corroboration.py`: Renders Figure 5d (`Visualisations/bubble_fea_final.pdf`) with concentric MCMC corroboration rings.

### 5. `flux_sampling/` — RHMC MCMC Uniform Flux Sampling (20 Paired CSMs)
- `scripts/`: RHMC uniform flux-sampling pipeline (`run_rhmc_new_chain.m`; polytope rounding followed by four independent chains per model [80 chains across the 20 models], Gelman-Rubin `R-hat` & ESS diagnostics, and Mann-Whitney U / Cliff's delta tests on `|v|`).
- `results/sampling_per_reaction.csv` & `sampling_summary.csv`: Reaction-level posterior flux statistics (`median |v|`, `mean |v|`, `R-hat`, `ESS`, `Cliff's delta`, BH `FDR`) across all 10 cancer-vs-normal pairs.
- `raw_chains/`: Local archive of raw thinned MCMC chain draws (git-ignored due to size).
