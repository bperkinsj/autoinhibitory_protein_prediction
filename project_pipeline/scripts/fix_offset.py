'''
This has to be its own step so that the offset can be fixed without re-downloading the files.
'''

import pandas as pd
import main

# Define the download path for the CIF files
cif_path = snakemake.input[0]
df_prot = pd.read_csv(snakemake.input[1], sep = '\t')

# Fix any offsets between the UniProt sequence and the PDB sequence in the CIF files
df_offsets = main.correct_offset(df_prot, cif_path)

df_offsets.to_csv(snakemake.output[0], sep = '\t', index = False)