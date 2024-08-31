# CSN Normalization and Comparison

This repository contains the implementation for reproducing the results of our paper:
"Comparison and Development of Cross-Study Normalization Methods for Inter-Species Transcriptional Analysis".

This project focuses on cross-study normalization methods for analyzing gene expression datasets from different species. We evaluate three existing normalization techniques—Empirical Bayes, Distance Weighted Discrimination, and Cross-Platform Normalization—on RNA-sequencing data. Additionally, we propose a novel normalization method designed to reduce experimental differences while preserving biologically significant variations. Our results demonstrate that our method offers improved conservation of biological differences compared to existing approaches, facilitating more accurate inter-species analysis.

## Table of Contents
1. [Data Preprocessing](#data-preprocessing)
2. [Additional Data](#additional-data)
3. [CSN Normalization](#csn-normalization)
4. [Other Normalization Methods](#other-normalization-methods)
5. [Estimation Method](#estimation-method)

## Data Preprocessing

Datasets were downloaded from GEO. The following experiments were used:
- GSE60424 (H1 - Human)
- GSE107011 (H2)
- GSE122597 (M1 - Mouse)
- GSE124829 (M2)

### Steps:
1. Download series matrices and place them in `data/raw_series`
2. Delete their first rows (file should start with the table itself)
3. Download the expression tables and place them in `data/raw_exp_tables`
4. Run `run_data_pre_process.m`
5. The data is now ready: 
   - Series matrices are in `data/series`
   - Expression tables are in `data/processed_exp_tables`

## Additional Data

### Info Files
The following files in `data/raw_series` are merges of two files downloaded from the experiment GEO site:
- `GSE60424_names.txt`
- `GSE107001_names.txt`
- `GSE124829_partition.xlsx`

Those files are merge of 2 files that were downloaded from the experiment GEO site:
GSExxxxx_series_matrix.txt.gz file and SraRunInfo.csv file that was downloaded from the SRA link in the experiment site


### Ensembl Data
- Orthologous gene tables were produced by Ensembl browser. 
- Original output can be found in `data/ortho_genes_tables`.
- Preprocessing code for these tables is in `prepare_ortho_tables_for_coordination.m`

Files `human_geneID_name.txt` and `mouse_geneID_name.txt` in `data/raw_exp_tables` were also produced using Ensembl.

## CSN Normalization

To run CSN normalization:
1. Place the preprocessed expression tables and series matrices in `CSN_normalization/data`.
2. Run `run_normalizationX.m` where X is 1 to 6.
3. Output will appear in `CSN_normalization/output`

## Other Normalization Methods

Code for other normalization methods mentioned in the paper is in the `other_norm_methods` directory.

Usage:
1. Place datasets in the same directory as the scripts.
2. Examples of how to run this code can be found in:
   - `runDWD.R`
   - `runEB.R`
   - `runXPNandBaselineScript.m`

## Estimation Method

The estimation method code is in the `estimation_method` directory.
The outputs of the normalization methods are currently stored in `estimation_method/data`.

###Steps:
1. Place series matrices, non-normalized and normalized datasets in `estimation_method/data`
2. Examples of running the estimation method can be found in:
   - `example_cross1.m`
   - `example_cross2.m`
   - `example_cross3.m`
   - `example_cross4.m`
   - `mouse_example.m`
   - `human_example.m`
