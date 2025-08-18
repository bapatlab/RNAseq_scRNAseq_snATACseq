# single-nucleus ATAC-seq

This folder contains the code and R environment to analyze the single-nucleus ATAC-seq data deposited in GSEXXX. Fragmet files were generated from raw fastq files using cellranger-atac count.

## Files
- `create_arrow_files_archR.ipynb`: create the arrow files for ArchR
- `basic_processing_archR.ipynb`: processing of ArchR files: DARs, UMAP, clustering, annotation, etc.
- `renv.lock`: R environment lock file for reproducibility.


## Description
- To reproduce the analysis, download the data from GEO and adjust the path in the notebook to the location of the data.
- Activate the R environment using `renv` and run the notebook.
