import sys
from pyprojroot.here import here
from pyhpo import Ontology, HPOSet, helper
import pandas as pd
import itertools
from scipy.spatial.distance import squareform
import h5py
import numpy as np

# Set basedir with here(), based on presence of .git file
here()
# Add script dir to path to import gcep.py and gcep_config.py
import os

print(here())

# Load HPO ontology object from HPO3
## Needs to occur before defining valid_hpo
Ontology()


def create_hpo_distance(hpo_set_list, kind = "omim", method = "graphic", combine = "BMA"):
    """
    Take a list of HPOSet objects and return a matrix of the distances between each set in the list

    Args:
        hpo_set_list (_type_): List of HPOSet objects objects 

    Returns:
        _type_: distance matrix
    """   
    
    # Create list of all tuples representing all pairwise combinations of HPO sets in HPO_set_list 
    hpoSet_combinations = [
        (a[0], a[1]) for a in itertools.combinations(hpo_set_list, 2)
    ]   

    # Get similarities of HPO set pairs from HPO3 
    # TODO parameterize kind, method, combine
    mtx_sim = helper.batch_set_similarity(
        hpoSet_combinations, 
        kind=kind,
        method=method,
        combine=combine
    )
    
    # Convert similarity matrix to distance matrix
    mtx_dist = squareform([1 - x for x in mtx_sim])
    return mtx_dist


def valid_hpo(hpo_str):
    """
    Simple function returning True/False based on whether and HPO
    string is in the HPO3 ontology

    Args:
        hpo_str (_type_): _description_

    Returns:
        _type_: _description_
    """    
    try:
        HPOSet.from_queries([hpo_str])
        return True
    except:
        return False
    
# Generate HPO table from query using hpoTable() method
df_probands = pd.read_csv()

# Filter to only valid HPOs
df_probands['valid'] = df_probands['HPO_ID'].map(valid_hpo) # Create column with validation status
df_probands = df_probands[df_probands['valid']].drop(columns=['valid']) # Filter to only valid HPOs and drop validation column
df_probands['label'] = [x.encode('ascii', 'ignore').decode('ascii') for x in df_probands['label'].to_list()]
print(df_probands)

# Get all unique HPO terms in df_probands and convert each HPO term to its own HPOset
## To be used to calculate distance between all individual HPO terms
all_hpos = [HPOSet.from_queries([x]) for x in set(df_probands['HPO_ID'].to_list())]

# Collapse all HPO terms from a given proband into a list
df_probands_long = df_probands
df_probands = df_probands.groupby(['Gene', 'Disease', 'label'])['HPO_ID'].apply(list).reset_index()
df_probands['HPO_ID'] = df_probands['HPO_ID'].map(HPOSet.from_queries)
# Create a concatenated ID column
df_probands['ID'] = df_probands['Gene'] + '__' + df_probands['Disease'] + '__' + df_probands['label']

# Create a dataframe of HPO metadata
df_hpo_meta = pd.DataFrame({
    'hpo_id': [x.terms()[0].id for x in all_hpos],
    'hpo_name': [x.terms()[0].name for x in all_hpos]
})

# Generate distance matrices for HPO terms
mtx_hpo_dist = create_hpo_distance(all_hpos, combine = combine)
# Generate distance matrix for probands
mtx_proband_dist = create_hpo_distance(df_probands['HPO_ID'].to_list(), combine = combine)

# Save all data to an HDF5 file
dt_str = datetime.now().strftime('%Y%m%d_%H%M%S')
hf_save_dir = here(f'data/clingen/api/{dt_str}')
os.makedirs(hf_save_dir, exist_ok=True)
date_range = f's{start_date.replace("-", "")}_e{end_date.replace("-", "")}'
hf_save_path = os.path.join(hf_save_dir, f'{active_gcep}_hpo_{date_range}_{combine}.h5')
print(f'Saving to {hf_save_path}')

with h5py.File(hf_save_path, 'w') as f:
    
    # Save HPO distance matrix and metadata
    ## HPO distance matrix
    f.create_dataset('hpo_distance', data=mtx_hpo_dist, compression = 'gzip')
    ## HPO metadata
    hpo_metadata_group = f.create_group('hpo_metadata')
    hpo_metadata_group.create_dataset('hpo_id', data=np.array(df_hpo_meta['hpo_id'].values, dtype='S'))
    hpo_metadata_group.create_dataset('hpo_name', data=np.array(df_hpo_meta['hpo_name'].values, dtype='S'))
    hpo_metadata_group.create_dataset('index', data=df_hpo_meta.index.values)
    hpo_metadata_group.attrs['columns'] = np.array(df_hpo_meta.columns.tolist(), dtype='S')
    
    # Save proband distance matrix and metadata
    ## Proband distance matrix
    f.create_dataset('proband_distance', data=mtx_proband_dist, compression = 'gzip')
    ## Proband metadata
    proband_metadata_group = f.create_group('proband_metadata')
    proband_metadata_group.create_dataset('gene', data=np.array(df_probands['Gene'].values, dtype='S'))
    proband_metadata_group.create_dataset('disease', data=np.array(df_probands['Disease'].values, dtype='S'))
    proband_metadata_group.create_dataset('proband_id', data=np.array(df_probands['label'].values, dtype='S'))
    proband_metadata_group.create_dataset('index', data=df_probands.index.values)
    proband_metadata_group.attrs['columns'] = np.array(df_probands.columns.tolist(), dtype='S')
    
    ## Original proband-hpo data
    proband_hpo_group = f.create_group('proband_hpo')
    proband_hpo_group.create_dataset('gene', data=np.array(df_probands_long['Gene'].values, dtype='S'))
    proband_hpo_group.create_dataset('disease', data=np.array(df_probands_long['Disease'].values, dtype='S'))
    proband_hpo_group.create_dataset('proband_id', data=np.array(df_probands_long['label'].values, dtype='S'))
    proband_hpo_group.create_dataset('hpo_id', data=np.array(df_probands_long['HPO_ID'].values, dtype='S'))
    proband_hpo_group.create_dataset('hpo_term', data=np.array(df_probands_long['HPO_term'].values, dtype='S'))
    proband_hpo_group.create_dataset('index', data=df_probands_long.index.values)
    proband_hpo_group.attrs['columns'] = np.array(df_probands_long.columns.tolist(), dtype='S')
    
    ## Gene metadata
    columns = ['Gene', 'Disease', 'Mode of Inheritance', 'Final Classification', 'Genetic Total Points', 'Experimental Total Points', 'Total Points']
    gene_metadata_group = f.create_group('gene_metadata')
    gene_metadata_group.create_dataset('gene', data=np.array(clingen_query.table['Gene'].values, dtype='S'))
    gene_metadata_group.create_dataset('disease', data=np.array(clingen_query.table['Disease'].values, dtype='S'))
    gene_metadata_group.create_dataset('inheritance', data=np.array(clingen_query.table['Mode of Inheritance'].values, dtype='S'))
    gene_metadata_group.create_dataset('classification', data=np.array(clingen_query.table['Final Classification'].values, dtype='S'))
    gene_metadata_group.create_dataset('points_genetic', data=np.array(clingen_query.table['Genetic Total Points'].values, dtype='S'))
    gene_metadata_group.create_dataset('points_experimental', data=np.array(clingen_query.table['Experimental Total Points'].values, dtype='S'))
    gene_metadata_group.create_dataset('points_total', data=np.array(clingen_query.table['Total Points'].values, dtype='S'))
    gene_metadata_group.create_dataset('index', data=clingen_query.table.index.values)
    gene_metadata_group.attrs['columns'] = np.array(columns, dtype='S')