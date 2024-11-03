import pickle

# Load the dictionaries from the pickle files
with open('chembl_protein_smile_results_dict.pkl', 'rb') as f:
    chembl_dict = pickle.load(f)

with open('old_protein_smile_results_dict.pkl', 'rb') as f:
    old_dict = pickle.load(f)

# Join the dictionaries
result_dict = {}
for key in chembl_dict:
    if key in old_dict:
        result_dict[key] = {**chembl_dict[key], **old_dict[key]}
    else:
        result_dict[key] = chembl_dict[key]

for key in old_dict:
    if key not in result_dict:
        result_dict[key] = old_dict[key]

print(result_dict)

# Save the result dictionary to a new pickle file
with open('protein_smile_results_dict_combined.pkl', 'wb') as f:
    pickle.dump(result_dict, f)