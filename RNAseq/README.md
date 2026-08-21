# bulk RNAseq data 

This folder contains the code and R environment to analyze the bulk RNAseq data deposited in PRJNA1212208 and shown in Extended Data Figures 7, 9 and 10. Count tables were generated using TEcount.

## Files
- `RNAseq_analysis.ipynb`: processing of RNAseq data
- `sample_metadata.csv`: metadata for the samples, showing sample_id, condition, group and replicates
- `renv.lock`: R environment lock file for reproducibility.


## Description
- To reproduce the analysis, download the data from GEO and adjust the path in the notebook to the location of the data.
- To restore the R environment from the lock file, run `renv::restore()` in the R terminal.
