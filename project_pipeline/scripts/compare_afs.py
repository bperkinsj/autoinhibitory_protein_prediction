'''
A script for getting the RMSD between given AlphaFold2 models. In the context of this pipeline, one model is the public model 
and the other is the model from the AlphaFold2 pipeline.
'''

import main
import csv
import pandas as pd

df = pd.read_csv(snakemake.input[0], sep='\t').astype('object')

# Drop unnecessary info
df = df[['uniprot', 'region_1', 'region_2', 'cluster', 'af_filename', 'cf_filename']].drop_duplicates().reset_index(drop=True)

path1 = snakemake.input[1] # Path to public models
path2 = snakemake.input[2] # Path to models from ColabFold pipeline
path3 = snakemake.input[3] # Path to save complex structures

rmsd_info = main.compare_af(df, path1, path2, path3)

with open(snakemake.output[0], 'w') as file:
    fields = ['uniprot', 'cluster', 'cf_filename', 'complex_rmsd', '1.0_aligned', '1.0_comp',
                '1.1_aligned', '1.1_comp', '1.2_aligned', '1.2_comp', '2.0_aligned', '2.0_comp',
                '2.1_aligned', '2.1_comp', '2.2_aligned', '2.2_comp', '2.3_aligned', '2.3_comp',
                '1_aligned', '1_comp', '2_aligned', '2_comp']
    writer = csv.DictWriter(file, fieldnames=fields, delimiter='\t')
    
    writer.writeheader()
    for item in rmsd_info:
        writer.writerow(item)