import numpy as np
import pandas as pd

from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity
from chembl_webresource_client.new_client import new_client


class Drug:
    def __init__(self, smiles: str, name=None, chembl_id=None, max_phase=None, indication=None, first_approval=None):
        # Get drug info from ChEMBL
        self.name = name
        self.chembl_id = chembl_id
        self.smiles = smiles
        self.max_phase = max_phase
        self.indication = indication
        self.first_approval = first_approval
        self.similarity = []
        self.gdsc_response = pd.DataFrame()

    def get_similar_drugs(self) -> list:
        similarity = new_client.similarity

        res = similarity.filter(smiles=self.smiles, similarity=40).only(['molecule_chembl_id', 'pref_name', 'max_phase', 'first_approval', 'indication_class', 'molecule_structures__canonical_smiles', 'similarity'])

        # Create similar drug objects and add ChEMBL info
        sim_drugs_list = [(Drug(drug["molecule_structures"]["canonical_smiles"],
                                name=drug["pref_name"],
                                chembl_id=drug["molecule_chembl_id"],
                                max_phase=drug["max_phase"],
                                indication=drug["indication_class"],
                                first_approval=drug["first_approval"]), drug["similarity"]) for drug in res]
        
        curr_drug = [drug for drug, sim in sim_drugs_list if (drug.smiles == self.smiles)][0]

        self.name = curr_drug.name
        self.chembl_id = curr_drug.chembl_id
        self.max_phase = curr_drug.max_phase
        self.indication = curr_drug.indication
        self.first_approval = curr_drug.first_approval
        
        self.similarity = [(drug["molecule_chembl_id"], drug["similarity"]) for drug in res]


        return sim_drugs_list
    
    def get_top_responders(self):
        gdsc_response = pd.read_excel("./data/GDSC2_fitted_dose_response_24Jul22.xlsx")
        name_variants = ["".join(self.name.split("-")), self.name.upper(), self.name.lower()]
        gdsc_response = gdsc_response[gdsc_response["DRUG_NAME"].isin(name_variants)]

        if len(gdsc_response) == 0:
            raise ValueError("Drug not found in GDSC2. Please select another.")
        
        self.gdsc_response = gdsc_response.sort_values("AUC")

    def smiles_similarity(mol1: str, mol2: str) -> float:
        # Create TfidfVectorizer object
        vectorizer = TfidfVectorizer()

        # Generate molecule vector
        mol1_tfidf = vectorizer.fit_transform([mol1])
        mol2_tfidf = vectorizer.fit_transform([mol2])

        # Calculate cosine similarity
        cosine_sim = cosine_similarity(mol1_tfidf, mol2_tfidf)

        return cosine_sim
    



