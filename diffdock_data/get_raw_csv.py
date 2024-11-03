import pandas as pd

# Load the original CSV file
df = pd.read_csv('test_smiles_single.csv')

# Create a new DataFrame by repeating each row 3 times with different protein names
protein_names = ['BRD4', 'sEH', 'HSA']
new_df = pd.concat([df.assign(protein_name=protein) for protein in protein_names])


# rename
# Save the new DataFrame to a new CSV file
new_df.to_csv('protein_smile_raw.csv', index=False)