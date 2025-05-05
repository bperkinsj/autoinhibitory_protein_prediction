from pdbecif.mmcif_io import CifFileReader, CifFileWriter
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import MMCIFParser, NeighborSearch, Selection
from pymol import cmd
import pymol
import pandas as pd
import numpy as np
from os.path import join
import requests
import os
import json

def make_dirs(*paths):
    for path in paths:
        try:
            os.mkdir(path)
        except FileExistsError:
            print('Folder already exists!')

def query_rcsb(uniprot_id, url):
    
    # Query test for pdb files associated with given UniProt accession number.
    query_text = {
    "query": {
      "type": "group",
      "logical_operator": "and",
      "nodes": [
        {
          "type": "terminal",
          "service": "text",
          "parameters": {
            "operator": "exact_match",
            "value": uniprot_id,
            "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession"
          }
        },
        {
          "type": "terminal",
          "service": "text",
          "parameters": {
            "operator": "exact_match",
            "value": "UniProt",
            "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name"
          }
        }
      ]
    },
    "request_options": {
      "return_all_hits": True
    },
    "return_type": "polymer_instance"
      }
      
    print("Querying RCSB PDB REST API for Uniprot_ID: %s" % uniprot_id)
      
    header = {'Content-Type': 'application/x-www-form-urlencoded'}
      
    response = requests.post(url, json=query_text)

    # In format 4CJ0 1P3I ...
    pdb_str = ''

    if response.status_code == 200:
          response_dic = response.json()
          for n in range(len(response_dic['result_set'])):
              pdb_str = pdb_str + response_dic['result_set'][n]['identifier'] + ' '
          
    else:
          pdb_str = np.nan
          print("Failed to retrieve results for Uniprot_ID: %s" % uniprot_id)

    return pdb_str


def prune_extra_chains(pdb_ids_str):

    #Turn the string into a list
    pdb_ids_w_chain = pdb_ids_str.strip().split(sep=' ')

    #Empty dictionary to fill.
    pdb_ids_dict = {}

    for pdb_id in pdb_ids_w_chain:

        #PDB ID (lowercase)
        pdb = pdb_id.split('.')[0].lower()

        #Chain label
        chain = pdb_id.split('.')[1]

        #Add the PDB ID as a key and the chain label as a value.
        if pdb not in pdb_ids_dict.keys():
            pdb_ids_dict[pdb] = [chain]
        else:
            pdb_ids_dict[pdb].append(chain)

    #Extract each PDB from the pdb_ids dict
    for pdb_id in pdb_ids_dict.copy():
            
        #Determine whether there are one or more chains.
        if len(pdb_ids_dict[pdb_id]) != 1:

            #If more than one chain, select the first chain as our representative.
            pdb_ids_dict[pdb_id] = pdb_ids_dict[pdb_id][0]

    # Now we convert them back to strings and add it all together.
    # Make a list of the values
    values_list = list(pdb_ids_dict.values())

    #Make a list of the keys
    key_list = list(pdb_ids_dict.keys())

    #Empty string to fill with my IDs.
    unique_pdb_ids = ''

    for n in range(len(pdb_ids_dict)):

        #Get the value
        chain = values_list[n][0]

        #Get the key
        key = key_list[n].lower()

        # Put them together
        pdb_id_chain_str = key + '.' + chain

         #Append to my unique_pdb_ids string
        unique_pdb_ids = unique_pdb_ids + ' ' + pdb_id_chain_str

    #Make the value of PDB at the index i equal to our new string.
    return unique_pdb_ids

def remove_chains(pdb_ids_chains):
    #The pdb ids will have their chains attached here (format example: 5ecy.A)
    pdb_ids_chains_list = pdb_ids_chains.split(sep=' ')

    #empty list to store pdb ids without chains
    pdb_ids_no_chains = []

    #Remove the chains from the PDB ids
    for pdb_id in pdb_ids_chains_list:
        pdb_id_only = pdb_id.split(sep='.')[0]
        pdb_ids_no_chains.append(pdb_id_only)
    
    return pdb_ids_no_chains

def expand_on_pdbs(df):
    # Convert PDB column to list-like
    df['pdb'] = df['pdb'].str.split(pat=' ')
    # Explode PDB column
    df = df.explode('pdb').reset_index(drop = True)
    # Split PDB ID and chain into separate columns
    df[['pdb', 'chain']] = df['pdb'].str.split(pat='.', expand = True)

    return df

def get_offset(fp, pdb, uniprot):
    # initiate reader object
    cfr = CifFileReader()
    cif_obj = cfr.read(fp, output='cif_wrapper')
    cif_data = list(cif_obj.values())[0]
    
    # Find the index of the chain of interest. If it doesn't exist, return 0.
    try:
        index = cif_data._struct_ref_seq.pdbx_db_accession.index(uniprot)

        # Extract the auth_seq start and db_seq start (from Uniprot) from the cif file
        pdb_start = int(cif_data._struct_ref_seq.seq_align_beg[index])
        unp_start = int(cif_data._struct_ref_seq.db_align_beg[index])
        offset = pdb_start - unp_start
    except ValueError:
        offset = 0

    print(f'Offset for {pdb}: {offset}')
    return offset

def fix_offset(pdb, fp, chain, offset):
    
    # Replace residue numbers in our chain of interest

    # Read in the CIF file
    cfr = CifFileReader()
    cif_obj = cfr.read(fp, output='cif_dictionary')
    # Convert to a Pandas DataFrame
    df = pd.DataFrame.from_dict(cif_obj[pdb.upper()]['_atom_site'])
    # Replace residue numbers
    for i in range(len(df)):
        if df.loc[i, 'label_asym_id']==chain:
            res_num = int(df.loc[i, 'label_seq_id'])
            new_res_num = res_num - offset
            df.loc[i, 'label_seq_id'] = str(new_res_num)
            df.loc[i, 'auth_seq_id'] = str(new_res_num)
        else:
            continue
    # Convert back to mmCIF-like dictionary
    cif_dict = df.to_dict(orient='list')
    cif_obj[pdb.upper()]['_atom_site'] = cif_dict

    # Write to file
    cfw = CifFileWriter(fp)
    cfw.write(cif_obj)

    return f'Successfully fixed {pdb}'
    
def string2range(x):
    
    """
    This function takes in a `string` representing a region of interest in a
    protein. The region of interest can be a single region or multiple regions
    of a protein. Returns a range for single regions or a list of ranges for
    multiple regions.
    
    Parameters:
    
        x (string): String containing a region or several regions of interest in a 
            protein.
            Format of x: single region -> 'start-end'
                         multiple regions -> 'start1-end1,start2-end2'
                     
    Returns:
    
        range or list of ranges: For single region proteins a range is returned. For 
            multiple region proteins a list of ranges is returned

            Format: single region -> range(start, end+1)
                    multiple region -> [range(start1, end1+1), range(start2, end2+1)]
    """
    # Handle instances with more than one range
    if ',' in x:
        list_temp = x.split(sep = ',') #list_temp = ['123-456,' '789-1111']
        for y in range(len(list_temp)): 
            list_temp[y] = list_temp[y].split(sep = '-') #list_temp[y] = [['123', '456'], ['789', '1111']]
        for y in range(len(list_temp)): 
            for x in range(len(list_temp[y])):
                list_temp[y][x] = int(list_temp[y][x]) #turns each list item into an integer

        # Make a range object with the bounds of the range. Note to the 
        # end a 1 has to be added in order to include the last position in the range
        for y in range(len(list_temp)): #[1, 2] where 1=[123, 456] and 2=[789, 1111]
            for x in range(len(list_temp[y])): #[123, 456]       
                list_temp[y] = list(range(list_temp[y][x], list_temp[y][x+1]+1)) #list_temp[0][0] = [123], list_temp[0][0+1]+1 or [456] + 1 = [457]
                break

        return list(set([item for sublist in list_temp for item in sublist]))

    # Handle instances with only one range
    else:
        list_temp = x.split(sep = '-')
        for y in range(len(list_temp)):
            list_temp[y] = int(list_temp[y]) #

        # Make a range object with the bounds of the region. Note to the 
        # end a 1 has to be added in order to include the last position in the range
        return list(range(list_temp[0], list_temp[1]+1))
    
def region_search_range(df):
    # Convert the domain region strings to ranges or lists of ranges
    df['region_1 search'] = df['region_1'].apply(lambda x: string2range(x))
    df['region_2 search'] = df['region_2'].apply(lambda x: string2range(x))

    return df
    
def get_structure_dict(name, fn, path):

    # Join the path and the file name
    full_path = join(path, fn)

    # To load a PDB file make a parser object
    parser = MMCIFParser(QUIET=True)
            
    # Then make a structure object
    structure = parser.get_structure(name, full_path)
            
    # Make an MMCIFDict object to grab more information form the .cif files
    mmcif_dict = MMCIF2Dict(full_path)

    return structure, mmcif_dict


def count_domain_residues(region1, region2, structure, label_chain):
    '''
    Count the number of residues in the domain and in the IAS.
    '''

     # Set all the counters to zero
    count_res = 0
    count_res_region_1 = 0
    count_res_region_2 = 0
    model_id = ''

    for model in structure:

        model_id = model.get_id()

        if model_id == 0:
        
            for chain in model:

                #Determine the current chain in the structure
                current_chain = chain.get_id()
                
                #Only act on the chain relevant to our protein of interest.
                if current_chain == label_chain:

                    print(f'We want {label_chain}. Currently analyzing {current_chain}.')
                        
                    # Get the model ID for later use
                    model_id = model.get_id()
                    # Get all the residues in the chain A
                    residues = chain.get_residues()
                    
                    # Iterate through all the residues in the chain and determine
                    # whether they belong to the IAS or to the Domain.
                    for residue in residues:
                        count_res = count_res + 1                    
                        
                        # Amino acid residues have an empty space in position zero
                        # of the id
                        if residue.get_id()[0] == ' ':
                            # The sequence position of the amino acid residue is stored
                            # in position 1 of the id
                            if residue.get_id()[1] in region1:
                                # print(residue.get_id()[1])
                                count_res_region_1 = count_res_region_1 + 1

                            elif residue.get_id()[1] in region2:
                                count_res_region_2 = count_res_region_2 + 1

            else:
                continue
        
        else:
            continue

    return count_res_region_1, count_res_region_2, count_res, model_id
            
def calculate_domain_completeness(region1, region2, count_in_region1, count_in_region2):
    # Calculate the percentage of residues in the IAS and in the Domain
    percent_in_region_1 = (count_in_region1/len(region1))*100 if len(region1) else 0
    
    percent_in_region_2 = (count_in_region2/len(region2))*100 if len(region2) else 0

    return percent_in_region_1, percent_in_region_2

def get_domain_residues(region1, region2, structure, label_model, label_chain):
    # Iterate through all the models in the structure
    for model in structure:
        
        # Analyze only the model that corresponds to the current row in df_prot
        if model.get_id() == label_model:
        
            for chain in model:
                
                # Analyze only the chain that corresponds to the current row in df_prot
                if chain.get_id() == label_chain:

                    print(f'Looking at chain {label_chain} for interface')
                    # Get all the atoms in the chain
                    atom_list = Selection.unfold_entities(chain, "A")
                    
                    # Make an empty list to store the atoms in the region_1 and in region_2
                    atoms_ns = []
                    
                    # Iterate through all the atoms in the chain and append the ones that occur inside
                    # the region_1 or the region_2 to the list
                    for atom in atom_list:
                        
                        # Get the parent residue for the atom
                        res = atom.get_parent()
                        
                        # Make sure the residue is an amino acid
                        if res.get_id()[0] == ' ':
                            
                            # Check whether the residue lies inside the region_1 or the region_2 and append
                            # the atoms of these residues into a list
                            if res.get_id()[1] in region1:
                                atoms_ns.append(atom)
                                
                            elif res.get_id()[1] in region2:
                                atoms_ns.append(atom)

                    return atoms_ns
                
def domain_neighborsearch(region1, region2, atoms):
    # Make an NeighborSearch object with all the atoms inside the region_1 and the region_2
    ns = NeighborSearch(atoms)
      
    # Search for all the interacting residues in the region_1 and in the region_2
    # with atoms that are within a 5 A radius 
    ns_all = ns.search_all(5, 'R')
    
    # Make a set to store the residues at the interface
    interface_res = set()
    
    # Save the interacting residue pairs as a list of tuples
    interacting_pairs = []
    
    # Iterate thorugh all the interacting residue pairs and append those that have a residue
    # in the region_1 and another in the region_2 to a list. Save the residue positions in a set
    for pairs in ns_all:
        
        res_0 = pairs[0].get_id()[1]
        res_1 = pairs[1].get_id()[1]
        
        if res_0 in region1 and res_1 in region2:
            interface_res.add(res_0)
            interface_res.add(res_1)
            interacting_pairs.append((res_0, res_1))
            
        elif res_1 in region1 and res_0 in region2:
            interface_res.add(res_0)
            interface_res.add(res_1)
            interacting_pairs.append((res_0, res_1))
            
    # Save the results in the appropriate columns of df_prot
    if len(interface_res) > 0 and len(interacting_pairs) > 0:
        return str(interacting_pairs), str(interface_res), len(interface_res)  
        
    else: 
        return np.nan, np.nan, np.nan
    
def to_string(x):
    x = str(x)
    
    if x == 'nan':
        return np.nan
    elif '{' in x:
        x = x.replace('{', '').replace('}', '').replace(' ', '')
        return x
    elif '[' in x:
        x = x.replace('[', '').replace(']', '').replace(' ', '')
        return x
    
def compare_atoms(gt_df, pred_df):

    present_atoms_pred = []
    extra_atoms_gt = []

    # Define atom_names for hydrogens
    hydrogens = ['HA', 'HB1', 'HB2', 'HB3', 'H', 'HA2', 'HA3', 'HG2', 'HG3', 'HD2', 'HD3', 'HE1', 'HE2',
                'HE3', 'HB', 'HG11', 'HG12', 'HG13', 'HG21', 'HG22', 'HG23', 'HE', 'HH11', 'HH12', 'HH21',
                'HH22', 'HE21', 'HE22', 'HD1', 'HZ', 'HH', 'HG1', 'HD21', 'HD22', 'HG', 'HD11', 'HD12', 
                'HD13', 'HD23', 'HZ1', 'HZ2', 'HZ3', 'HH2']
    # Define possible alternate conformations
    alt_locations = ['B', 'C', 'D', 'E']

    for atom in range(len(gt_df)):
        # Define rows to be skipped (hydrogens, alternate conformations, or extra NMR models)
        if gt_df.loc[atom, 'label_atom_id'] in hydrogens or gt_df.loc[atom, 'label_alt_id'] in alt_locations or gt_df.loc[atom, 'pdbx_PDB_model_num'] != '1':
            extra_atoms_gt.append(atom)
            continue

        # Define minimum  parameters to select unique rows
        gt_atom_name = gt_df.loc[atom, 'label_atom_id']
        gt_residue_name = gt_df.loc[atom, 'label_comp_id']
        gt_residue_number = gt_df.loc[atom, 'label_seq_id']

        pred_row = pred_df.loc[(pred_df['label_atom_id'] == gt_atom_name) 
                               & (pred_df['label_seq_id'] == gt_residue_number) 
                               & (pred_df['label_comp_id'] == gt_residue_name)]
        if pred_row.empty != True:
            present_atoms_pred.append(pred_row.index)
        else:
            extra_atoms_gt.append(atom)

    return present_atoms_pred, extra_atoms_gt

def drop_unshared_atoms(gt_df, pred_df, present_atoms_pred, extra_atoms_gt):
    # Select all rows in pred not present in gt
    total_atoms = list(pred_df.index)
    na_atoms_array = np.setdiff1d(total_atoms, present_atoms_pred)
    na_atoms = sorted(na_atoms_array)

    # Create new pred data frame exclusively with atoms present in gt
    pred_trim = pred_df.drop(index=na_atoms)
    gt_trim = gt_df.drop(index=extra_atoms_gt)

    return gt_trim, pred_trim

def assert_equal_size(gt_trim_df, pred_trim_df):
    try:
        assert len(pred_trim_df) == len(gt_trim_df)
        return True
    except AssertionError:
        gt_sim = gt_trim_df.drop(['id', 'Cartn_x', 'Cartn_y', 'Cartn_z', 'occupancy', 'B_iso_or_equiv', 'type_symbol', 'pdbx_formal_charge'], axis=1)
        pred_sim = pred_trim_df.drop(['id', 'Cartn_x', 'Cartn_y', 'Cartn_z', 'occupancy', 'B_iso_or_equiv', 'type_symbol', 'pdbx_formal_charge'], axis=1)
        diff = pd.concat([gt_sim, pred_sim]).drop_duplicates(keep=False)
        diff.to_csv('./data/AssertionError.tsv', sep='\t')
        print(diff)
        print('AssertionError! Check file')
        return False
    
def trim_stats(uniprot, pdb, gt, gt_trim, pred, pred_trim):
    gt_perc = len(gt_trim) / len(gt)
    pred_perc = len(pred_trim) / len(pred)
    trim_values_dict = {'uniprot': uniprot,
                        'pdb': pdb,
                        'gt_len': len(gt),
                        'gt_trim_len': len(gt_trim),
                        'pred_len': len(pred),
                        'pred_trim_len': len(pred_trim),
                        'gt_perc': gt_perc,
                        'trim_perc': round(pred_perc, 3)}
    
    return trim_values_dict
    
def create_region_dict(region, region_num):
    '''Create a dictionary containing an ID for every region in the domain
    For instance, if domain 1 has 123-222,333-444, then make dict {1.0: 123-222+333-444, 1.1: 123-222, 1.2: 333-444}.
    #.0 always contains the full number of regions.'''
    full_region = region.strip()
    # Each dict will have at least one entry.
    region_dict = {f'{region_num}.0': full_region}
    if ',' in region:
        full_region = region.replace(',', '+')
        # Substitute full_region in the dict to one with + in place of ,
        region_dict[f'{region_num}.0'] = full_region
        regions_list = region.split(',')
        for i in range(len(regions_list)):
            # Make subregions 1-indexed to preserve #.0 as full.
            subregion = i + 1
            region_dict[f'{region_num}.{subregion}'] = regions_list[i]

    return region_dict
    
def load_and_select(gt_fn, pred_fn, region_1, region_2):
    # Load and select native and pred pdbs
    cmd.delete('all')
    cmd.load(gt_fn, 'native')
    cmd.load(pred_fn, 'pred')

    for obj in ['native','pred']:
        # select region1 and region2
        for key in region_1:
            # example: native_1.1, native and resi 111-222
            resi_range = region_1[key]
            cmd.select(f'{obj}_{key}', f'{obj} and resi {resi_range}')
        for key in region_2:
            resi_range = region_2[key]
            cmd.select(f'{obj}_{key}', f'{obj} and resi {resi_range}')

def align_and_calculate(align_reg_key, comp_region_key):
    # superimpose aligned region and calculate two rmsds: One for aligned region and one for complementary region (for example, rmsd for "aligned" region 1.1 
    # and complementary region 2.0)
    rmsds = []
    try:
        align = cmd.align(f'native_{align_reg_key}', f'pred_{align_reg_key}')
        rmsd = cmd.rms_cur(f'native_{align_reg_key}', f'pred_{align_reg_key}', matchmaker=4)
        rmsds.append(round(rmsd, 3))
        rmsd = cmd.rms_cur(f'native_{comp_region_key}', f'pred_{comp_region_key}', matchmaker=4)
        rmsds.append(round(rmsd, 3))
        return rmsds

    except pymol.CmdException:
        print(f'Region {align_reg_key} missing')
        rmsds = [-1, -1]
        return rmsds
    
def get_region_averages(rmsds):
    '''
    If there are multiple regions in a domain, calculate their average rmsd.
    '''

    for item in rmsds:
        if item['1.0_aligned'] == 0:
            if item['1.1_aligned'] == -1:
                item['1_aligned'] = item['1.2_aligned']
                item['1_comp'] = item['1.2_comp']
            elif item['1.2_aligned'] == -1:
                item['1_aligned'] = item['1.1_aligned']
                item['1_comp'] = item['1.1_comp']
            else:
                item['1_aligned'] = round((item['1.1_aligned'] + item['1.2_aligned']) / 2, 3)
                item['1_comp'] = round((item['1.1_comp'] + item['1.2_comp']) / 2, 3)
        else:
            item['1_aligned'] = item['1.0_aligned']
            item['1_comp'] = item['1.0_comp']

    for item in rmsds:
        if item['2.0_aligned'] == 0 and item['2.3_aligned'] == 0:
            if item['2.1_aligned'] == -1:
                item['2_aligned'] = item['2.2_aligned']
                item['2_comp'] = item['2.2_comp']
            elif item['2.2_aligned'] == -1:
                item['2_aligned'] = item['2.1_aligned']
                item['2_comp'] = item['2.1_comp']
            else:
                item['2_aligned'] = round((item['2.1_aligned'] + item['2.2_aligned']) / 2, 3)
                item['2_comp'] = round((item['2.1_comp'] + item['2.2_comp']) / 2, 3)
        elif item['2.0_aligned'] == 0 and item['2.3_aligned'] != 0:
            if item['2.1_aligned'] == -1:
                item['2_aligned'] = round((item['2.2_aligned'] + item['2.3_aligned']) / 2, 3)
                item['2_comp'] = round((item['2.2_comp'] + item['2.3_comp']) / 2, 3)
            elif item['2.2_aligned'] == -1:
                item['2_aligned'] = round((item['2.1_aligned'] + item['2.3_aligned']) / 2, 3)
                item['2_comp'] = round((item['2.1_comp'] + item['2.3_comp']) / 2, 3)
            elif item['2.3_aligned'] == -1:
                item['2_aligned'] = round((item['2.1_aligned'] + item['2.2_aligned']) / 2, 3)
                item['2_comp'] = round((item['2.1_comp'] + item['2.2_comp']) / 2, 3)
            else:
                item['2_aligned'] = round((item['2.1_aligned'] + item['2.2_aligned'] + item['2.3_aligned']) / 3, 3)
                item['2_comp'] = round((item['2.1_comp'] + item['2.2_comp'] + item['2.3_comp']) / 3, 3)
        else:
            item['2_aligned'] = item['2.0_aligned']
            item['2_comp'] = item['2.0_comp']

    return rmsds

def common_member(a, b):
    a_set = set(a)
    b_set = set(b)

    if (a_set & b_set):
        return list((a_set & b_set))
    else:
        return []
    
def region_bounds(x):
    """
    This function takes in a `string` representing a region of interest in a
    protein. The region of interest can be a single region or multiple regions
    of a protein. Returns the bounds of the regions in list form.
    
    Parameters:
    
        x (string): String containing a region or several regions of interest in a 
            protein.
            Format of x: single region -> 'start-end'
                         multiple regions -> 'start1-end1,start2-end2'
                     
    Returns:
    
        region boundaries in list form

            Format: single region -> [start, end+1]
                    multiple region -> [[start1, end1+1], [start2, end2+1]]
    """
    # Handle instances with more than one range
    if ',' in x:
        list_temp = x.split(sep = ',') #list_temp = ['123-456,' '789-1111']
        for y in range(len(list_temp)): 
            list_temp[y] = list_temp[y].split(sep = '-') #list_temp[y] = [['123', '456'], ['789', '1111']]
        for y in range(len(list_temp)): 
            for x in range(len(list_temp[y])):
                list_temp[y][x] = int(list_temp[y][x]) #turns each list item into an integer

        # Make a range object with the bounds of the range. Note to the 
        # end a 1 has to be added in order to include the last position in the range
        for y in range(len(list_temp)): #[1, 2] where 1=[123, 456] and 2=[789, 1111]
            for x in range(len(list_temp[y])): #[123, 456]       
                list_temp[y] = [list_temp[y][x], list_temp[y][x+1]+1] #list_temp[0][0] = [123], list_temp[0][0+1]+1 or [456] + 1 = [457]
                break

        return list_temp

    # Handle instances with only one range
    else:
        list_temp = x.split(sep = '-')
        for y in range(len(list_temp)):
            list_temp[y] = int(list_temp[y]) #

        # Make a range object with the bounds of the region. Note to the 
        # end a 1 has to be added in order to include the last position in the range
        return [[list_temp[0], list_temp[1]+1]]
    
def pae_from_json(path, fn):
    '''Read in the json file, which is in the format:
    [{"predicted_aligned_error":[[0, 1, 3, 5, 19, ...], [0, 4, 12, 38, ...], ...]}]
    '''

    try:
        f = open(os.path.join(path, fn))
        data = json.load(f)
        pae = data['predicted_aligned_error']
        array = np.array(pae)
    except FileNotFoundError:
        print(f'File {fn} not found')
        array = np.empty((1,1,))
        array[:] = np.nan

    return array


def calculate_pae_mean(prot_array, reg_a, reg_b):
    '''
    Gives the mean pae for all regions of interest compared against all regions of interest (reg1 to reg1, reg1 to reg2, reg2 to reg2)
    Reg_a and Reg_b are given as arrays.
    '''

    '''
    Method proceeds like this:
    Let's say we're handed reg1 from a protein, which is "1-2, 3-4". It's given to us in the form:
    reg_a = [[1, 2], [3, 4]], reg_b = [[1, 2], [3, 4]]]
    So we perform:
    mean(prot_array[1:3, 1:3]) + mean(prot_array[1:3, 3:5]) + mean(prot_array[3:5, 1:3]) + mean(prot_array[3:5, 3:5]). Then take mean of all that.
    '''

    # First is reg1 on reg1
    means = []
    for i in range(len(reg_a)):
        for n in range(len(reg_b)):
            a_start = reg_a[i][0]
            a_end = reg_a[i][1]
            b_start = reg_b[n][0]
            b_end = reg_b[n][1]
            sub_array = prot_array[a_start:a_end, b_start:b_end]
            sub_mean = np.mean(sub_array)
            means.append(sub_mean)

    mean = np.mean(means)

    return (mean)

def uniprot_dirs(*path_list, uniprot):
    "Make subdirectories for the given uniprot id"
    for path in path_list:
        try:
            os.mkdir(path + '/' + uniprot + '/')
        except:
            pass

def cif_to_pdb(*fns):
    # Change a file from .cif to .pdb and return the new file name

    for f in fns:
    #Change gt file to pdb
        cmd.load(f, 'protein')
        cmd.save(f.replace('.cif', '.pdb'), selection='protein')
        cmd.delete('all')

    return f.replace('.cif', '.pdb')

def add_AF_filename(df, fp):
    # Specifically for AlphaFold files because CF have clusters

    # Get list of filenames from the file path
    filenames = os.listdir(fp)

    # Keep only the .cif files
    filenames = [f for f in filenames if f.endswith('.cif') or ('rank_001' in f and f.endswith('.pdb'))]

    # Determine which files are in the dataframe
    df_proteins = df['uniprot'].tolist()

    # Create a dictionary of the filenames and the corresponding uniprot
    filename_dict = {}
    for filename in filenames:
        if filename.endswith('.cif'):
            uniprot = filename.split('-')[1]
            filename_dict[uniprot] = filename
        elif filename.endswith('.pdb'):
            uniprot = filename.split('_')[0]
            filename_dict[uniprot] = filename

    # Create a list of the filenames in the order of the dataframe
    AF_filenames = []
    for protein in df_proteins:
        AF_filenames.append(filename_dict[protein])

    # Add the list to the dataframe
    df['af_filename'] = AF_filenames

    return df

def select_regions(gt_fn, pred_fn, region_1, region_2):
    cmd.delete('all')
    cmd.load(gt_fn, 'native')
    cmd.load(pred_fn, 'pred')

    # Region format is currently "123-456,789-1111". We need to
    # change it to "123-456+789-1111" to select both regions at once.
    region_1 = region_1.replace(',', '+')
    region_2 = region_2.replace(',', '+')

    # Make a dict to iterate through
    regions = {'1': region_1, '2': region_2}

    for obj in ['native','pred']:
        for key in regions:
            resi_range = regions[key]
            cmd.select(f'{obj}_{key}', f'{obj} and resi {resi_range}')

def alter_chain(gt_fn, pred_fn):

    # First set the native and pred chains to C
    cmd.alter('native', 'chain="C"')
    cmd.alter('pred', 'chain="C"')

    # Then change region 1 to chain B
    cmd.alter('native_1', 'chain="B"')
    cmd.alter('pred_1', 'chain="B"')

    # And change region 2 to chain A
    cmd.alter('native_2', 'chain="A"')
    cmd.alter('pred_2', 'chain="A"')

    # Save the files
    cmd.save(gt_fn, 'native')
    cmd.save(pred_fn, 'pred')

def add_CF_filename(df, path):

    dirs = os.listdir(path)

    file_dict = {'uniprot': [], 'cluster': [], 'cf_filename': []}

    # Create a new dataframe with the file names 
    for d in dirs:
        files = os.listdir(os.path.join(path, d))
        for f in files:
            if '.pdb' in f:
                uniprot = f.split('_')[0] # filename example: P28482_U10-000_unrelaxed_rank_001_alphafold2_multimer_v2_model_1_seed_000.pdb
                cluster = f.split('_')[1]
                rank = f.split('_')[4]
                if rank != '001': # We want the best structure
                    continue

                # Populate the dictionary
                file_dict['uniprot'].append(uniprot)
                file_dict['cluster'].append(cluster)
                file_dict['cf_filename'].append(f)

    # Create a new dataframe from the dictionary
    file_df = pd.DataFrame(file_dict)

    # Merge the new dataframe with the old one
    fn_df = pd.merge(df, file_df, on=['uniprot'], how='left')

    fn_df = fn_df.dropna(subset=['cf_filename'])

    return fn_df

def get_pdb_struct_dict(name, fn, path):
    # Join the path and the file name
    full_path = os.path.join(path, fn)

    # To load a PDB file make a parser object
    parser = PDBParser(QUIET=True)
            
    # Then make a structure object
    structure = parser.get_structure(name, full_path)

    return structure

def adjust_trimmed_bounds(reg1, reg2):
    '''
    This function takes in two numpy arrays of region boundaries in a sequence 
    (either [1, 4] or [[1, 5], [8, 12], ...]), finds the minimum from both
    lists, and subtracts the minimum boundary from every value to adjust for 
    trimmed sequences.'''

    min_bound = min(reg1.min(), reg2.min())

    result1 = np.subtract(reg1, min_bound - 1)
    result2 = np.subtract(reg2, min_bound - 1)

    return result1, result2