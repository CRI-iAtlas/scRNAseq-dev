# scRNAseq-dev: Code to prepare single-cell RNA-seq datasets for ingress in iAtlas


## Step 1: Create a Synapse folder for storing original and processed data

To ensure reproducibility of the code and improve the preparation of manifestts for the iAtlas database, the original data and processed files are stored in Synapse.

The datasets listed in this repo have their associated data stored at: https://www.synapse.org/Synapse:syn52077487/files/ (you need permission to access this folder)

We recommend the creation of a new folder with a title that references the dataset, and the following subfolders:
- Original files: files from the original source. This includes direct exports from h5ad files (for example, a csv with UMAP coordinates extracted from a h5ad object).
- iAtlas files: files created during data processing.


## Step 2: format a h5ad file (data_prep.ipynb files)

If the original data is available in an h5ad file, 

### Requirements
- adata.obs.cell_type_iatlas 
- adata.obs.sample_id
- adata.obs.patient_id
- normalized counts, computed with sc.pp.log1p(adata)

We also export the following data into a csv file in Synapse:
- UMAP coordinates
- OBS table


