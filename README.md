# Autoinhibitory Protein Prediction
This repository contains the main component of the pipeline for comparing AlphaFold2-generated predictions to experimental structures.

## Contents

### environments
Contains the necessary environments to run the pipeline and Jupyter notebooks. Requires [Anaconda](https://www.anaconda.com/download).

### project_pipeline
The pipeline used for analyzing high-throughput structural data can be found here. Please see below for instructions on how to run it. Within project_pipeline/data/starting_proteins/ are the lists of initial proteins (autoinhibitory, two-domain or multi-domain, and obligate) analyzed in the paper. ```fifteen_ai_proteins_clusters.csv``` contains the experimental files and AF-Cluster-generated clusters compared per protein.

### source_protein_lists
Contains the original autoinhibitory protein lists provided by Dr. Dokyun Na and supplemented by Jorge A. Holguin-Cruz and Brooks Perkins-Jechow.

### multi_domain.ipynb
Contains code to collect the list of multi-domain proteins from the PDB.

### paper_figs.ipynb
Contains code to generate the bulk of the paper figures.

### paper_lineage_score.ipynb
Contains code to generate the Shannon entropy and lineage score figures.

### structural_clustering.ipynb
Contains code to perform hierarchical clustering of the experimental structures and generate the plot.

## The pipeline

### Install the environment
Make sure you have [Anaconda](https://www.anaconda.com/download) installed.

To install the ```rmsd_snek``` environment, follow the documentation on [creating an environment from an environment.yml file](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file).

### Running the pipeline
The first step of the pipeline is the ```Snakefile``` in the ```project_pipeline``` folder. The pipeline requires a tab-separated (tsv) file in the format

| uniprot | region_1 | region_2 |
| ---     | ---      | ---      | 
| P28583  | 15-139   | 211-412  | 

where "uniprot" is the UniProt ID, "region_1" is the sequence range of the inhibitory module, and "region_2" is the sequence range of the functional domain. It also requires all of the AlphaFold-predicted files to be placed in ```data/Alphafold_cif```. Once those are provided, the pipeline can be run as per normal Snakemake.

First, navigate to the project_pipeline folder:

$ cd ./project_pipeline/

and then run Snakemake with the following command:

$ snakemake -s Snakefile --cores

If you need to adjust the filename for the file you're passing to Snakemake, you can go into the Snakefile (```./project_pipeline/Snakefile```) and under rule pdb_ids change the input file 'data/proteins.tsv' to whatever your file is (i.e. 'data/proteins_to_measure.tsv') or such. Make sure that the file is tab-separated.

The final output will be to the file ```./project_pipeline/data/rmsds.tsv```.

## Re-creating figures

