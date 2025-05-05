'''This script takes in the files with autoinhibited and two-domain protein clusters
and calculates the mean predicted aligned error for each cluster.'''

import pandas as pd
import main

ai = pd.read_csv('./data/ai_pdb_clusters.tsv', sep='\t')
md = pd.read_csv('./data/md_pdb_clusters.tsv', sep='\t')

ai_path = 'data/input/Colabfold_cif/autoinhibited/'
md_path = 'data/input/Colabfold_cif/multi_domain/'

affix = ''
suffix = 'predicted_aligned_error_v1.json'

# Get rid of some redundancy because we only need the uniprots and clusters
ai = ai[['uniprot', 'cluster', 'region_1', 'region_2']].drop_duplicates().reset_index(drop=True) 
md = md[['uniprot', 'cluster', 'region_1', 'region_2']].drop_duplicates().reset_index(drop=True)

# Calculate mean PAE for each cluster
ai_pae = main.mean_paes(ai, ai_path, affix, suffix, cluster=True)
md_pae = main.mean_paes(md, md_path, affix, suffix, cluster=True)

# Save the files
ai_pae.to_csv('./data/ai_cluster_pae.tsv', sep='\t', index=False)
md_pae.to_csv('./data/md_cluster_pae.tsv', sep='\t', index=False)