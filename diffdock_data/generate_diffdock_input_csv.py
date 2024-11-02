import pandas as pd
import random
import sys

# read K from command line
if len(sys.argv) > 1:
    try:
        K = int(sys.argv[1])
    except ValueError:
        K = None
else:
    K = 10

# Load the CSV file
df = pd.read_csv('belka/test.csv')

if K is not None:
    # Select K random ids
    ids = random.sample(range(1, len(df)), K)
    bindings = df.iloc[ids]
else:
    # read the ids from ids.txt
    with open('ids.txt', 'r') as f:
        ids = f.read().splitlines()

    # convert the ids to integers
    ids = [int(i) for i in ids]
    bindings = df[df['id'].isin(ids)]

bindings = bindings[['id', 'molecule_smiles', 'protein_name']]

# change binding column names to match the expected column names
bindings.columns = ['complex_name', 'ligand_description', 'protein_path']

# convert protein_name to protein_path
bindings['protein_path'] = '../diffdock_data/proteins/' + bindings['protein_path'] + '.pdb'

# convert id to string and prepend 'complex_'
bindings['complex_name'] = 'complex_' + bindings['complex_name'].astype(str)

# Add empty columns for protein_sequence
bindings['protein_sequence'] = ''

# write the selected lines to a new CSV file
bindings.to_csv('diffdock_input.csv', index=False)
