'''
A script to add cluster information and filenames to our dataframes.
'''

import os
import utils
import pandas as pd

'''
Retrieve clusters and filenames.
'''
df = pd.read_csv(snakemake.input[0], sep='\t')

path = snakemake.input[1]

# Drop rmsd info from full-depth to experimental comparisons, plus other unnecessary info
df = df[['uniprot', 'pdb', 'region_1', 'region_2', 'af_filename', 'chain']]

# Add the filenames
df_clusters = utils.add_CF_filename(df, path)

# Save the file
df_clusters.to_csv(snakemake.output[0], sep='\t')