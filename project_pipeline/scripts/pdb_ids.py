# -*- coding: utf-8 -*-
"""
Created on Thu May 14 11:18:27 2020

@author: Jorge Holguin

Edited by Brooks Perkins-Jechow
"""

import pandas as pd
import utils
import main

# Define the download path for the CIF files
pdb_path = snakemake.input[0]
af_path = snakemake.input[1]
df_prot = pd.read_csv(snakemake.input[2], sep = '\t')
# df_prot = pd.read_csv('../data/protein_list.tsv', sep = '\t')

print('Querying RCSB for PDB IDs.')

df_prot = main.get_pdb_ids(df_prot)
print('Successfully retrieved IDs. Proceeding to download structures.')

# Download the pdb files
main.download_pdb_files(df_prot, pdb_path)

# Make a new dataframe with each PDB ID in a separate row and chains in their own column
df_pdb = utils.expand_on_pdbs(df_prot)

# Add filenames
df_pdb = utils.add_AF_filename(df_pdb, af_path)
df_pdb['gt_fn'] = df_pdb['pdb'] + '.cif'

# Save the dataframe as a tsv file
df_pdb.to_csv(snakemake.output[0], sep = '\t', index = False)

