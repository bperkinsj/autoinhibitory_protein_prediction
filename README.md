# Autoinhibitory Protein Prediction

[![DOI](https://zenodo.org/badge/978359796.svg)](https://doi.org/10.5281/zenodo.17307831)


## Abstract
Many proteins operate by toggling between distinct conformations, yet most structure predictors remain fixated on a single static fold. We benchmarked AlphaFold2, AlphaFold3, and recent variants on autoinhibited proteins, a type of allosterically regulated protein whose sequences encode at least two functional states. AlphaFold consistently assigned lower confidence to these flexible proteins, indicating it has internalized aspects of energy landscapes. AlphaFold3 showed the highest accuracy when supplied with post-translational modifications or binding partners, correctly modeling both functional states in ~50\% of cases—though this relied on prior knowledge.  Without such cues, uniform MSA subsampling, not signal deconvolution, was unexpectedly more effective in capturing alternative conformations, particularly when using diverse sequences. We also found that predicted aligned error (PAE) outperformed pLDDT in assessing model quality. Overall, we propose the use of full-depth and uniformly subsampled MSAs combined with low PAE selection for known autoinhibited proteins,  which yields higher-confidence predictions, often favoring compact, closed conformations that may offer insight into regulatory mechanisms. By analyzing over 100 previously untested proteins, our study also highlights key challenges in predicting structures shaped by complex, evolved energy landscapes.

## Contents
This repository contains the main component of the pipeline for comparing AlphaFold2-generated predictions to experimental structures.

### environments
Contains the necessary environments to run the pipeline and Jupyter notebooks. Requires [Anaconda](https://www.anaconda.com/download).

### project_pipeline
The pipeline used for analyzing high-throughput structural data can be found here. Please see below for instructions on how to run it. Within project_pipeline/data/starting_proteins/ are the lists of initial proteins (autoinhibitory, two-domain or multi-domain, and obligate) analyzed in the paper. ```fifteen_ai_proteins_clusters.csv``` contains the experimental files and AF-Cluster-generated clusters compared per protein.

### source_protein_lists
Contains the original autoinhibitory protein lists provided by Dr. Dokyun Na and supplemented by Jorge A. Holguin-Cruz and Brooks Perkins-Jechow.

| file | description |
| ---  | ---         |
| Autoinhibited proteins (more added by Brooks) | Original curated list provided by Dr. Dokyun Na, with additional proteins manually added by Brooks Perkins-Jechow (see paper for methods). Contains annotations of evidence for autoinhbiition in selected proteins. |
| Jorge_Algonquin_Autoinhibition_Proteins_w_Structures | Curated list provided by Jorge Holguin-Cruz. Provides references to sourced literature. |
| autoinhibited_proteins_no_autoinhibited_structures_dokyun_na | List of autoinhibited proteins with no known autoinhibited structures in the PDB. Provided by Dr. Dokyun Na. |
| autoinhibited_proteins_no_structures_dokyun_na | List of autoinhibited proteins with no known structures in the PDB. Provided by Dr. Dokyun Na. |

### multi_domain.ipynb
Contains code to collect the list of multi-domain proteins from the PDB.

### paper_figs.ipynb
Contains code to generate the bulk of the paper figures.

### paper_lineage_score.ipynb
Contains code to generate the Shannon entropy and lineage score figures.

### structural_clustering.ipynb
Contains code to perform hierarchical clustering of the experimental structures and generate the plot.

## The Pipeline

### Part 1: Install the environment
- Make sure you have [Anaconda](https://www.anaconda.com/download) installed.

- To install the ```rmsd_snek``` environment, follow the documentation on [creating an environment from an environment.yml file](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file).
  * Please note the [dependencies](https://github.com/bperkinsj/autoinhibitory_protein_prediction/blob/main/environments/rmsd_snek.yml) that are required.


### Part 2: Running the pipeline

1. Check file format.

Verify that your input file is in the correct format for the ```Snakefile``` in the ```project_pipeline``` folder. The pipeline requires a tab-separated (tsv) file in the format shown below:

| uniprot | region_1 | region_2 |
| ---     | ---      | ---      | 
| P28583  | 15-139   | 211-412  | 

where "uniprot" is the UniProt ID, "region_1" is the sequence range of the inhibitory module, and "region_2" is the sequence range of the functional domain. It also requires all of the AlphaFold-predicted files to be placed in ```data/Alphafold_cif```. Once those are provided, the pipeline can be run as per normal Snakemake.

2. Navigate to the project pipeline folder

$ cd ./project_pipeline/

3. Run Snakemake

$ snakemake -s Snakefile --cores

The final output will be to the file ```./project_pipeline/data/rmsds.tsv```.

#### Additional notes

If you need to adjust the filename for the file you're passing to Snakemake, you can go into the Snakefile (```./project_pipeline/Snakefile```) and under rule pdb_ids change the input file 'data/proteins.tsv' to whatever your file is (i.e. 'data/proteins_to_measure.tsv') or such. Make sure that the file is tab-separated.

## Re-creating figures

All plots can be re-created by running all cells in paper_figs.ipynb, paper_lineage_score.ipynb, and structural_clustering.ipynb. Composite figures were created in Inkscape alongside some manual editing of axis labels for a few figures such as the number of base or alternate conformations.

## Issues
If you have any questions or problems please leave a comment under the Issues tab.
