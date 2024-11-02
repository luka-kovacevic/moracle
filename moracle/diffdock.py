from typing import Tuple
import pickle

def get_diffdock(protein_name, smile) -> Tuple[str, float]:
    '''
    Given a protein name and a smile string, return the location of a zip file
    containing the results of DiffDock analysis, i.e. 
        1. The protein structure as a PDB file
        2. The ligand structure as a SDF file
    '''
    
    # Load the precomputed dictionary
    with open('../diffdock_data/protein_smile_results_dict.pkl', 'rb') as f:
        protein_smile_dict = pickle.load(f)

    # Check if the protein name is in the dictionary
    if protein_name not in protein_smile_dict:
        return None
    
    if smile not in protein_smile_dict[protein_name]:
        return None
    
    experiment_id, confidence = protein_smile_dict[protein_name][smile]
    
    return f'../diffdock_data/zips/{experiment_id}.zip', confidence
    
if __name__ == '__main__':
    protein_name = 'BRD4'
    smile = 'Cc1onc(-c2ccccc2)c1C(=O)N1CCC[C@@H](c2cccc(-c3ccc(C4CCCC4)cc3)c2)[C@@H]1C(=O)N[Dy]'
    print(get_diffdock(protein_name, smile))