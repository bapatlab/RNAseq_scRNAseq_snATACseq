# bulk RNAseq data 

This folder contains the code and R environment to analyze the bulk RNAseq data deposited in PRJNA1212208. Count tables were generated using TEcount.

## Files
- `LFDvsHFD_Ear_bulk_RNA-seq.R`: analysis to create heatmap shown in Figure 3a
- `bulkRNAseq.ipynb`: analysis to convert Count Tables to a count matrix & create volcano plot (Extended Data Figure 9), heatmap and PCA plot (Extended Data Figure 12)
- `renv.lock`: R environment lock file for reproducibility.


## Description
- To reproduce the analysis, download the data from GEO and adjust the path in the notebook to the location of the data.
- To restore the R environment from the lock file, run `renv::restore()` in the R terminal.
