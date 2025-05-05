'''
A script for calculating average pLDDT scores of the entire protein and all annotated domains for given proteins.
'''

import pandas as pd
import main

# Define file paths
df_ai_fd = pd.read_csv(snakemake.input[0], sep='\t').astype('object')
df_md_fd = pd.read_csv(snakemake.input[1], sep='\t').astype('object')
df_sd_fd = pd.read_csv(snakemake.input[2], sep='\t').astype('object')
df_ai_cl = pd.read_csv(snakemake.input[3], sep='\t').astype('object')
df_md_cl = pd.read_csv(snakemake.input[4], sep='\t').astype('object')  # TODO: Pass mean_plddt one dataframe containing only uniprots and af_filenames,   # and another (unchanged) dataframe containing cf_filenames. 
ai_fd_fp = snakemake.input[5]
md_fd_fp = snakemake.input[6]
sd_fd_fp = snakemake.input[7]
ai_c_fp = snakemake.input[8]
md_c_fp = snakemake.input[9]

# Calculate average pLDDT scores for AlphaFold2 models
ai_fd_scores = main.mean_plddt(df_ai_fd, ai_fd_fp)
md_fd_scores = main.mean_plddt(df_md_fd, md_fd_fp)
sd_fd_scores = main.mean_plddt_single_domain(df_sd_fd, sd_fd_fp)

# Calculate average pLDDT scores for ColabFold models
ai_c_scores = main.mean_plddt(df_ai_cl, ai_c_fp, fnt='cf_filename')
md_c_scores = main.mean_plddt(df_md_cl, md_c_fp, fnt='cf_filename')

# Write results to file
ai_fd_scores.to_csv(snakemake.output[0], sep='\t', index=False)
md_fd_scores.to_csv(snakemake.output[1], sep='\t', index=False)
sd_fd_scores.to_csv(snakemake.output[2], sep='\t', index=False)
ai_c_scores.to_csv(snakemake.output[3], sep='\t', index=False)
md_c_scores.to_csv(snakemake.output[4], sep='\t', index=False)