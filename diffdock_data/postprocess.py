import pandas as pd
import pickle
import zipfile
import os

# Paths
csv_path = 'protein_smile_raw.csv'
pkl_path = 'protein_smile_results_dict.pkl'
proteins_dir = 'proteins/'
results_dir = 'results/'
zips_dir = 'zips/'

# Read the protein_smile_results_dict.pkl
with open(pkl_path, 'rb') as pkl_file:
    protein_smile_results_dict = pickle.load(pkl_file)

# Read the protein_name, molecule_smiles combinations from the CSV using pandas
df = pd.read_csv(csv_path)

for index, row in df.iterrows():
    protein_name = row['protein_name']
    molecule_smiles = row['molecule_smiles']

    # Get the corresponding id from the dictionary
    id, confidence = protein_smile_results_dict.get(protein_name, {}).get(molecule_smiles)
    if id is None:
        continue
    id = str(id)

    # Paths to the protein and molecule files
    protein_file = os.path.join(proteins_dir, f'{protein_name}.pdb')

    # search for the molecule file in the results directory it is of the form <id>/rank1_confidence{some number}.sdf
    # we want to find the value of {some_number}
    for file in os.listdir(f'{results_dir}complex_{id}'):
        if file.endswith('.sdf') and 'rank1_confidence' in file:
            confidence = file.split('confidence')[1].split('.sdf')[0]
            break
    
    #access that file
    molecule_file = os.path.join(f'{results_dir}complex_{id}', f'rank1_confidence{confidence}.sdf')

    # update the confidence value in the dictionary

    if confidence is not None:
        confidence = float(confidence)
    else:
        continue
    protein_smile_results_dict[protein_name][molecule_smiles] = (id, confidence)

    # Create a zip file containing the protein and molecule files
    zip_file_path = os.path.join(zips_dir, f'{id}.zip')
    with zipfile.ZipFile(zip_file_path, 'w') as zipf:
        zipf.write(protein_file, os.path.basename(protein_file))
        zipf.write(molecule_file, os.path.basename(molecule_file))

# Write the updated dictionary back to the pickle file
with open(pkl_path, 'wb') as pkl_file:
    pickle.dump(protein_smile_results_dict, pkl_file)