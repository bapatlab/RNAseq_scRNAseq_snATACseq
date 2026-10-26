# single-nucleus ATAC-seq

This folder contains the code and R environment to analyze the single-nucleus ATAC-seq data (Figure 4i, j, k and Extended Data Figure 10). Fragmet files were generated from raw fastq files using cellranger-atac count.

## Files
- `01_create_arrows.ipynb`: create the arrow files for ArchR
- `02_analysis.ipynb`: processing of ArchR files: DARs, UMAP, clustering, annotation, etc.
- `renv.lock`: R environment lock file for reproducibility.
- Files inside `inputs` folder: HFD and LFD barcodes.
- Files inside `refernece` folder: mm39 annotations, blacklist and genes, as well as checksums for Arrow files


## Description
- To reproduce the analysis, download the data from GEO and adjust the path in the notebook to the location of the data.
- Activate the R environment using `renv` and run the notebook.
