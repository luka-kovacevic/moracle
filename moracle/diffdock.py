from typing import Tuple
import pickle
import pandas as pd

predictions_csv = '../data/predictions.csv'

def get_binding_prob(protein_name, smile) -> float:
    '''
    Given a protein name and a smile string, return the binding probability
    (as caluclated by the model from Kaggle)
    '''
    
    # Load the predictions csv with pandas
    df = pd.read_csv(predictions_csv)

    # find the row with the protein name and smile
    row = df[(df['protein_name'] == protein_name) & (df['molecule_smiles'] == smile)]

    # if the row does not exist or two rows exist, return None
    if len(row) != 1:
        return None
    
    # return the binding probability
    prob = row['binds'].values[0]

    return prob

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
        return None, None
    
    if smile not in protein_smile_dict[protein_name]:
        return None, None
    
    experiment_id, confidence = protein_smile_dict[protein_name][smile]
    
    return f'../diffdock_data/zips/{experiment_id}.zip', confidence
    
if __name__ == '__main__':
    protein_name = 'BRD4'
    smile = 'C#CCCC[C@H](Nc1nc(Nc2ccc(F)c(OC)c2)nc(Nc2ccc(F)c(OC)c2)n1)C(=O)N[Dy]'
    print(get_diffdock(protein_name, smile))
    print(get_binding_prob(protein_name, smile))

    smile = 'COc1cc(C#N)ccc1-c1ccc2c(c1)c(C(=O)N[Dy])cn2[C@H]1CCCN(C(=O)C2CCCC2(F)F)C1'
    print(get_diffdock(protein_name, smile))
    print(get_binding_prob(protein_name, smile))
