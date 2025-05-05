'''
This script fixes several PDB files that had the incorrect UniProt ID or
db_align_beg. Although not part of the pipeline, this script was used in
the generation of the original data, and so is necessary to run to reproduce the results.
'''


from pdbecif.mmcif_io import CifFileReader
from pdbecif.mmcif_tools import MMCIF2Dict
from pdbecif.mmcif_io import CifFileWriter, CifFileReader
import pandas as pd
import numpy as np
from Bio.PDB import MMCIFParser, NeighborSearch, Selection
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

df = pd.read_csv('./data/file_fixes.csv', sep='\t').astype('object')


for i in range(len(df)):
    uniprot = df.loc[i, 'uniprot']
    pdb = df.loc[i, 'pdb']
    fp = './data/input/RCSB_cif/' + uniprot + '/' + pdb + '.cif'
    seq = df.loc[i, 'seq']

    # Read in the file
    cfr = CifFileReader()
    cif_obj = cfr.read(fp, output='cif_dictionary')
    struct = cif_obj[pdb.upper()]['_struct_ref_seq']
    
    if len(struct['align_id']) > 1:
        if not pd.isnull(df.loc[i, 'old_uniprot']):
            old_uniprot = df.loc[i, 'old_uniprot']
            # Fix old uniprot to new uniprot:
            index = struct['pdbx_db_accession'].index(old_uniprot)
            struct['pdbx_db_accession'][index] = uniprot

                # Correct db_align_beg
            if not pd.isnull(df.loc[i, 'seq']):
                struct['db_align_beg'][index] = str(int(seq))

        else:
            index = struct['pdbx_db_accession'].index(uniprot)
            struct['db_align_beg'][index] = str(int(seq))
    else:
        if not pd.isnull(df.loc[i, 'old_uniprot']):
            old_uniprot = df.loc[i, 'old_uniprot']
            # Fix old uniprot to new uniprot:
            struct['pdbx_db_accession'] = uniprot

            if not pd.isnull(df.loc[i, 'seq']):
                struct['db_align_beg'] = str(int(seq))

        else:
            struct['db_align_beg'] = str(int(seq))

    # Add back to df_cif
    cif_obj[pdb.upper()]['_struct_ref_seq'] = struct

    # Write to file
    cfw = CifFileWriter(fp)
    cfw.write(cif_obj)

    print(f'Successfully fixed {pdb}')