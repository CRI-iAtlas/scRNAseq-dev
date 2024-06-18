# scRNAseq-dev: Code to prepare single-cell RNA-seq datasets for ingress in iAtlas

scRNAseq-dev stores code used to prepare new single-cell RNA-seq (scRNAseq) datasets for further processing for intake in the CRI iAtlas database.

This is the **first stage** to add a new scRNAseq dataset in iAtlas. Currently, there are **4 stages** in this process, summarized below: 
1. Data preparation (code in this repo), with the following outputs:

    a. h5ad file with normalized expression as a layer, and harmonized cell type names

    b. UMAP coordinates in a csv file

2. Compute gene expression statistics *(currently in notebooks stored in this repo, to be converted to a new workflow)*, with the following output:

    a. A csv file with statistics of frequency and average of expression of a given gene accross cells in a cell type

3. Run workflow to compute pseudobulk and geneset scoring (https://github.com/CRI-iAtlas/nf-pseudobulk), with the following output:

    a. Pseudobulk values of expression for each combination of sample + cell type

    b. GSEA results using values from (3a) and genesets described in [Thorsson et al, 2018](https://gdc.cancer.gov/about-data/publications/panimmune).

4. Create manifest files for upload in the CRI iAtlas database (https://github.com/CRI-iAtlas/iatlasManifests)


In this stage, we recommend the steps described below to prepare the data to be compatible with downstream processing.

## Step 1: Create a Synapse folder for storing original and processed data

To ensure reproducibility of the code and improve the preparation of manifests for the iAtlas database, the original data and processed files are stored in Synapse.

The datasets listed in this repo have their associated data stored at: https://www.synapse.org/Synapse:syn52077487/files/ (you need permission to access this folder). You can use a different folders or use another strategy to store these files, but changes to loading and uploading data in the files will be necessary.

We recommend the creation of a new folder with a title that references the dataset, and the following subfolders:
- Original files: files from the original source. This includes direct exports from h5ad files (for example, a csv with UMAP coordinates extracted from a h5ad object).
- iAtlas files: files created during data processing.


## Step 2: format a h5ad file (data_prep.ipynb files)

If the original data is available in an h5ad file, 

### Requirements
- adata.obs.cell_type_iatlas 
- adata.obs.sample_id
- adata.obs.patient_id
- gene symbol
- normalized counts, computed with sc.pp.log1p(adata)

We also export the following data into a csv file in Synapse:
- UMAP coordinates
- OBS table


