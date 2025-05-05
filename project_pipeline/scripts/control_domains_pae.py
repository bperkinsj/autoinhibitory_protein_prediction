'''
This script calculates the mean predicted aligned error for the subset of proteins labeled as "single domain" by
SCOP-e on the PDB that have only one domain annotated on UniProt and the subset of proteins labeled as "multi domain" by
SCOP-e on the PDB that have two domains annotated on Uniprot. 
See autoinhibition_protein_data/single_domain.ipynb for how these proteins were retrieved.
'''

import pandas as pd
import main

# Read in the file for the single domain proteins
df = pd.read_csv('./data/single_domain_domains.csv', sep=',').astype('object')
# Define file path for Alphafold predicted aligned error files
af_path = './data/input/Alphafold_single_domain/'
affix = 'AF-'
suffix = '-F1-predicted_aligned_error_v4.json'

# Calculate mean PAE for each protein
df_mean_pae = main.mean_pae_single_domain(df, af_path)

# Save file
df_mean_pae.to_csv('./data/single_domain_pae.tsv', sep='\t', index=False)

# Read in the file for the multi domain proteins
df = pd.read_csv('./data/multi_domain.tsv', sep='\t').astype('object')

# Define file path for Alphafold predicted aligned error files
af_path = './data/input/Alphafold_multi_domain/'

# Calculate mean PAE for each protein
df_mean_pae_md = main.mean_paes(df, af_path, affix, suffix)

df_mean_pae_md = df_mean_pae_md.dropna(subset='mean_pae_1_1')

# Save file
df_mean_pae_md.to_csv('./data/multi_domain_pae.tsv', sep='\t', index=False)
