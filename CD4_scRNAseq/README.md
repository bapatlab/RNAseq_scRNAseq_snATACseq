# CD4+ T cell single-cell RNA sequencing

This folder contains the code and python environment to analyze the single-cell RNA sequencing data of sorted CD4+ T cells deposited in GSEXXX and shown in Figure 2, Extended Data Figure 7, and Extended Data Figure 8.

## Files
- `cd4_sorted_analysis.ipynb`: QC, normalization, clustering, annotation, etc.
- `scrnaseq_environment.yml`: python environment to run the code
- `Combined_Gene_List.xlsx`: list of genes used to calculate Th17 pathogenicity score


## Description
- To reproduce the analysis, download the data from GEO and adjust the path in the notebook to the location of the data.
- Activate the conda environment (scrnaseq_environment.yml)
- Run the notebook
