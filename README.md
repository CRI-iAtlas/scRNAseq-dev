# scRNAseq-dev: Code to prepare single-cell RNA-seq datasets for ingress in iAtlas

scRNAseq-dev stores code used to prepare new single-cell RNA-seq (scRNAseq) datasets for further processing for intake in the CRI iAtlas database.

This is the **1st stage** to add a new scRNAseq dataset in iAtlas. Currently, there are **4 stages** in this process, summarized below: 
1. Data preparation (code in this repo), with the following outputs:

    a. h5ad file with normalized expression as a layer, and harmonized cell type names

    b. UMAP coordinates in a csv file

2. Compute gene expression statistics *(currently in notebooks stored in this repo, to be converted to a new workflow)*, with the following output:

    a. A csv file with statistics of frequency and average of expression of a given gene accross cells in a cell type

3. Run workflow to compute pseudobulk and geneset scoring (https://github.com/CRI-iAtlas/nf-pseudobulk), with the following output:

    a. Pseudobulk values of expression for each combination of sample + cell type

    b. GSEA results using values from (3a) and genesets described in [Thorsson et al, 2018](https://gdc.cancer.gov/about-data/publications/panimmune).

4. Create manifest files for upload in the CRI iAtlas database (https://github.com/CRI-iAtlas/iatlasManifests). Please note that this stage involves clinical mapping between the original data source and the structure in iAtlas, and this process is still heavily manual.


In this 1st stage, we recommend the steps described below to prepare the data to be compatible with downstream processing.

## Step 1: Create a Synapse folder for storing original and processed data

To ensure reproducibility of the code and improve the preparation of manifests for the iAtlas database, the original data and processed files are stored in Synapse.

The datasets listed in this repo have their associated data stored at: https://www.synapse.org/Synapse:syn52077487/files/ (you need permission to access this folder). You can use a different folder or another strategy to store these files, but changes to loading and uploading data in the files will be necessary.

We recommend the creation of a new folder with a title that references the dataset, and the following subfolders:
- Original files: files from the original source. This includes direct exports from h5ad files (for example, a csv with UMAP coordinates extracted from a h5ad object).
- iAtlas files: files created during data processing.


## Step 2: inspect the original data (data_prep.ipynb notebooks)

The majority of the scRNAseq studies currently included in iAtlas have its data in a h5ad (AnnData) file. For those that provide data in a different format, we recommend organizing it in this file format, as the code in stages 2 and 3 (pseudobulk, GSEA and cell expression stats) use functions specific for h5ad.

For a dataset to be eligible for inclusion in iAtlas, we currently require the following data:
- summarized cell x gene level expression values
- clinical annotation with: cell id, sample id, patient id, and cell annotation

After this data is gathered, we create a Jupyter Notebook exploring the following topics:
1. Check clinical annotation, as in several cases it includes information that will be necessary in the clinical mapping to iAtlas.
2. Check the different levels of cell type annotation. We use the annotattion set by cell ontology.
3. Check the QC metrics. At this moment, we are using data that was already processed by its authors, so we need to make sure that the data we are working with is the one that was already filtered. 

## Step 3: format the h5ad file (data_prep.ipynb notebooks)

The following should be present to the h5ad file for it to be ready for the next stages of processing.

### Requirements
- adata.obs.cell_type_iatlas: in case the cell type annotation is not compatible with the one used in iAtlas, create the `cell_type_iatlas` annotation in the obs table to add the harmonized version.
- expression matrix using gene symbols.
- normalized counts, computed with `scanpy.pp.log1p(adata)` and stored as a layer called `normalized`.

We recommend uploading this h5ad file in the Synapse folder created in Step 1, under the folder "iAtlas files".

## Step 4: export UMAP and obs data as csv (data_prep.ipynb notebooks)

When available, we also export the following data into a csv file in Synapse:
- UMAP coordinates
- OBS table
These tables are useful when creating the manifest files (Stage 4), and having them as csv files speed up the process.

We recommend uploading it in the Synapse folder created in Step 1, under the folder "iAtlas files".