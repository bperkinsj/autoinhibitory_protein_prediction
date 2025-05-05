from biopandas.mmcif import PandasMmcif
import pandas as pd

df = pd.read_csv('./data/proteins_pdb_best.tsv', sep='\t').astype('object')

entity_ids = []
for i in range(len(df)):

    # select PDB, file path, and author chain
    pdb = df.loc[i, 'PDB ID']
    cif_path = (f'./data/input/RCSB_cif_best/{pdb}.cif')
    auth_chain = df.loc[i, 'Auth_chain']

    # Initialize PandasMmcif object and load mmcif file
    pmcif = PandasMmcif()
    _ = pmcif.read_mmcif(cif_path)
    pred = pmcif.df['ATOM']

    #Filter to only our author chain and select entity id
    entity_id_df = pred.loc[pred['auth_asym_id']==auth_chain, 'label_entity_id']
    entity_id = entity_id_df.iloc[0]
    entity_ids.append(entity_id)

df.insert(12, 'Entity_id', entity_ids)

df.to_csv('./data/proteins_pdb_best.tsv', sep='\t', index=False)