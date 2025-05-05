'''
This is a script for taking the trimmed PDB files and splitting them into two chains
for use in DockQ. The autoinhibitory region will become chain B, while the functional domain
and the rest of the protein will remain chain A.
'''
import pandas as pd
import utils
import main

df = pd.read_csv(snakemake.input[0], sep='\t')
gt_path_in = snakemake.input[1]
pred_path_in = snakemake.input[2]
gt_path_out = snakemake.output[0]
pred_path_out = snakemake.output[1]

# If the source files are from Colabfold clusters, then we need "cluster=True" in main.split_chains
result = 'Colabfold_pdb_trim' in pred_path_in
alphafold = 'Alphafold' in gt_path_in
# Make directories
utils.make_dirs(gt_path_in, pred_path_in, gt_path_out, pred_path_out)

# If we are comparing the full-depth structures to the cluster structures, then we need to add the full-depth filenames (af_filename)
if alphafold:
    uniprots = df[['uniprot']].drop_duplicates().reset_index(drop=True)

    af_fn = utils.add_AF_filename(uniprots, gt_path_in)

    df = pd.merge(df, af_fn, on='uniprot', how='left')

# Split chains
main.split_chains(df, gt_path_in, pred_path_in, gt_path_out, pred_path_out, cluster=result, pred_only=alphafold)