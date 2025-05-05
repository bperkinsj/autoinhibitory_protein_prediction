'''
Determine the interface residues of the ColabFold-generated structures
'''
import pandas as pd
import main

path = snakemake.input[1]

df = pd.read_csv(snakemake.input[0], sep='\t')

# Drop unnecessary info
df = df[['uniprot', 'cluster', 'region_1', 'region_2', 'cf_filename']].drop_duplicates().reset_index(drop=True)

# Get the interfaces of the proteins
print('Finding interfaces...')
interacting = main.get_af_interfaces(df, path, cluster=True)

# Save the dataframe
interacting.to_csv(snakemake.output[0], sep='\t', index=False)