import pandas as pd
import main
import utils
import csv

gt_in_path = snakemake.input[0]
gt_trim_path = snakemake.input[1]
pred_in_path = snakemake.input[2]
pred_trim_path = snakemake.input[3]


# Read in the reference dataframe
df_prot = pd.read_csv(snakemake.input[4], sep = '\t')

# Make the directories
utils.make_dirs(gt_in_path, gt_trim_path, pred_in_path, pred_trim_path)

# Trim the files to make sure any mutated residues are not accidentally counted in our domain completeness
trim_values, df_prot = main.trim_cf_pdb(df_prot, gt_in_path, gt_trim_path, pred_in_path, pred_trim_path,
                                        gt_format = '{uniprot}/{pdb}.cif', pred_format = '{uniprot}/{filename}')

# Save the trim values
df_prot.to_csv(snakemake.output[0], sep = '\t', index = False)