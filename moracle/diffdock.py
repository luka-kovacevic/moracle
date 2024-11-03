from typing import Tuple
import pickle
import pandas as pd
import random

predictions_csv = '../data/predictions.csv'
# predictions_csv = '../data/predictions_dummy.csv'


def get_binding_prob(protein_name, smile) -> float:
    '''
    Given a protein name and a smile string, return the binding probability
    (as caluclated by the model from Kaggle)
    '''
    
    # Load the predictions csv with pandas
    df = pd.read_csv(predictions_csv)

    # find the row with the protein name and smile
    row = df[(df['protein_name'] == protein_name) & (df['molecule_smiles'] == smile)]

    # if the row does not exist or two rows exist, return a random number between 0 and 0.14
    if len(row) != 1:
        return random.uniform(0, 0.14)
    
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
    with open('../diffdock_data/protein_smile_results_dict_combined.pkl', 'rb') as f:
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

    # malformed smile so definitely won't be in the dictionary
    smile = '[[[[[COc1cc(C#N)ccc1-c1ccc2c(c1)c(C(=O)N[Dy])cn2[C@H]1CCCN(C(=O)C2CCCC2(F)F)C1'
    print(get_diffdock(protein_name, smile))
    print(get_binding_prob(protein_name, smile))


    #check demo zips
    demo_smiles = {'molecule1': 'C#CCCC[C@H](Nc1nc(Nc2ccc(F)c(OC)c2)nc(Nc2ccc(F)c(OC)c2)n1)C(=O)N[Dy]',
                   'molecule1mod': 'C#CCCC[C@@H](C(N[Dy])=O)Nc1nc(Nc2cc(OO)c(O)cc2)nc(Nc2cc(OC)c(F)cc2)n1',
                   'molecule4': 'C#CCCC[C@H](Nc1nc(NCc2ccc(C)cc2N2CCCC2)nc(Nc2ccc(F)c(OC)c2)n1)C(=O)N[Dy]'}
    for smile_name in ['molecule1', 'molecule1mod', 'molecule4']:
        for protein_name in ['BRD4', 'sEH', 'HSA']:
            smile = demo_smiles[smile_name]

            print(f'{smile_name}_{protein_name}, {get_diffdock(protein_name, smile)}, {get_binding_prob(protein_name, smile)}')

    #check demo2 zips
    demo_smiles = {'chembl1': 'COc1cc(Nc2nc(CO)cc(N(C)C)n2)ccc1F',
                   'chembl2': 'COc1cc(Nc2nc(C(C)C)cc(N(C)C)n2)ccc1F',
                   'chembl3': 'COc1cc(Nc2nc(C)cc(N(C)C)n2)ccc1F'}

    for smile_name in ['chembl1', 'chembl2', 'chembl3']:
        for protein_name in ['BRD4', 'sEH', 'HSA']:
            smile = demo_smiles[smile_name]

            print(f'{smile_name}_{protein_name}, {get_diffdock(protein_name, smile)}, {get_binding_prob(protein_name, smile)}')