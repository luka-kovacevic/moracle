import pandas as pd
import random
import sys
import pickle

# read K from command line
if len(sys.argv) > 1:
    try:
        K = int(sys.argv[1])
    except ValueError:
        K = None if sys.argv[1] != 'demo' else 'demo'
else:
    K = 10

# Load the CSV file
df = pd.read_csv('belka/test.csv')

# # load the precomputed dictionary
# with open('protein_smile_results_dict.pkl', 'rb') as f:
#     protein_smile_dict = pickle.load(f)

if K is None:
    # Select K random ids
    ids = random.sample(range(1, len(df)), K)
    protein_smile = df.iloc[ids]
    protein_smile = protein_smile[['molecule_smiles', 'protein_name', 'id']]

else:
    # read in the protein-smile combinations from protein_smile_raw.csv
    protein_smile = pd.read_csv('protein_smile_raw.csv')

    # add an id column to the dataframe
    if K == 'demo':
        demo_ids = []
        for protein in ['BRD4', 'sEH', 'HSA']:
            for name in ['molecule1', 'molecule1mod', 'molecule4']:
                demo_ids.append(f'{name}_{protein}')
            
        protein_smile['id'] = demo_ids
    else:
        protein_smile['id'] = range(1, len(protein_smile) + 1)


# populate protein_smile_results_dict with ids, None
protein_smile_results_dict = {}
for i in range(len(protein_smile)):
    protein_name = protein_smile['protein_name'].to_list()[i]
    smile = protein_smile['molecule_smiles'].to_list()[i]
    id = protein_smile['id'].to_list()[i]
    if protein_name not in protein_smile_results_dict:
        protein_smile_results_dict[protein_name] = {}
    protein_smile_results_dict[protein_name][smile] = (id, None)

# write the protein_smile_results_dict to a pickle file
with open('protein_smile_results_dict.pkl', 'wb') as f:
    pickle.dump(protein_smile_results_dict, f)


# change binding column names to match the expected column names
protein_smile.columns = ['ligand_description', 'protein_path', 'complex_name']

# convert protein_name to protein_path
protein_smile['protein_path'] = '../diffdock_data/proteins/' + protein_smile['protein_path'] + '.pdb'

# convert id to string and prepend 'complex_'
protein_smile['complex_name'] = 'complex_' + protein_smile['complex_name'].astype(str)

# Add empty columns for protein_sequence
protein_smile['protein_sequence'] = ''

# write the selected lines to a new CSV file
protein_smile.to_csv('diffdock_input.csv', index=False)
