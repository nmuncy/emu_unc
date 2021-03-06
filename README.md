# emu_unc

Scripts and data used to investigate the interaction of clinical diagnosis and tract with memory performance, amygdaloid activity during emotional judgments, and amygdaloid connectivity (PPI) with nodes in emotional regulation systems.

Seven main steps occurred when processing and analyzing data:
1. Pre-processing functional data via [https://github.com/emu-project/func_processing](https://github.com/emu-project/func_processing).
2. Additional pre-processing of functional data and deconvolution via `cli/func1_tfs.sh`, deconvolution workflow from step (1).
3. ROI construction and coefficient extraction via `cli/func2_roi.sh`. 
4. Functional PPI analysis via `cli/func3_ppi.sh` and `cli/func4_ppi_roi.sh`.
5. Diffusion pre-processing as reported [here](https://www.sciencedirect.com/science/article/pii/S221315822200002X).
6. Tract modeling with AFQ via `cli/diff1_run_afq.sh`.
7. Statistical tests and anlayses via `stats/manuscripts_stats.R`.

A description of repo directories below.

## cli
Wrapper scripts for various scripts found in `resources`. Numbered for respective order in their pipelines.

## data
Raw data generated by the workstream, and formatted dataframes used in analyses.

* `tract_profiles.csv` - raw output generated by AFQ.
* `AFQ_dataframe.csv` - formatted dataframe of AFQ, with extra participant information included.
* `Coef_*` - raw ROI and PPI coefficients extracted during Study scene valence judgments.
* `df_*_` - formatted dataframes of ROI, PPI coefficients. 

### timing_files  
Directory containing timing and reference JSON files needed by [https://github.com/emu-project/func_processing](https://github.com/emu-project/func_processing) (`workflow.control_afni.control_deconvolution`).
  
## env
A few documents for environment configuration.

## resources
Support functions (R) and work scripts (shell, python) for diffusion, functional, and statistical analyses. Also contains original pipeline used before restructuring of the project.

## stats
R scripts for making dataframes from `data/df_*`, conducting stats, and simulating data for examples.
