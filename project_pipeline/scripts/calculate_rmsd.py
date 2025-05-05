import pandas as pd
import csv
import main
import utils

gt_in_path = snakemake.input[0]
pred_in_path = snakemake.input[1]
complex_path = snakemake.input[2]
df = pd.read_csv(snakemake.input[3], sep='\t').astype('object')

# Make the directories
utils.make_dirs(gt_in_path, pred_in_path, complex_path)

# Calculate and save the rmsd info
rmsd_info = main.get_rmsds(df, gt_in_path, pred_in_path, complex_path)

with open(snakemake.output[0], 'w') as file:
    fields = ['uniprot', 'pdb', 'complex_rmsd', '1.0_aligned', '1.0_comp',
                '1.1_aligned', '1.1_comp', '1.2_aligned', '1.2_comp', '2.0_aligned', '2.0_comp',
                '2.1_aligned', '2.1_comp', '2.2_aligned', '2.2_comp', '2.3_aligned', '2.3_comp',
                'region_1', 'region_2', 'percent_region_1', 'percent_region_2',
                'chain', 'gt_fn', 'af_filename', 'complex_fn', '1_aligned', '1_comp',
                '2_aligned', '2_comp']
    writer = csv.DictWriter(file, fieldnames=fields, delimiter='\t')
    
    writer.writeheader()
    for item in rmsd_info:
        writer.writerow(item)