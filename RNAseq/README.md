# bulk RNAseq data 

This folder contains the code and R environment to analyze the bulk RNAseq data deposited in GSEXXX. Count tables were generated using TEcount.

## Files
- `TE_analysis.ipynb`: analysis of TEcount output, annotates genes, creates an overview of the percentage of each RNA class and analyzes retrovirus composition and frequency.
- `TE_DESeq2.ipynb`: basic bulk RNAseq analysis, including DESeq2, heatmaps, volcano plots, etc. 
- `renv.lock`: R environment lock file for reproducibility.


## Description
- To reproduce the analysis, download the data from GEO and adjust the path in the notebook to the location of the data.
- Activate the R environment using `renv` and run the notebook.
