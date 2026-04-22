## Protocol for identification of DNA Binding Factor enrichment in chromatin accessibility data to define a Persister Cell Signature in High-Grade Serous Ovarian Carcinoma 

The pcs_protocol.R file contains all the code necessary to run the analysis.


### Required data:
1. Counts matrix with peak data from a snATAC-seq or paired multiome data in h5 format. You can also start directly with a Seurat/Signac object with an ATAC assay.
2. Grouping metadata for the cells in the same format as in sample_data/metadata.csv
3. ReMap2022 non-redundant peaks BED file, and optional addons.bed with additional peaks in the same format as sample_data/addons.bed

