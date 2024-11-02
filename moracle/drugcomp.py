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

    def get_similar_drugs(self, clin_phase=4) -> list:
        similarity = new_client.similarity
        res = similarity.filter(smiles=self.smiles, similarity=40).filter(max_phase=clin_phase).only(['molecule_chembl_id', 'pref_name', 'max_phase', 'first_approval', 'indication_class', 'molecule_structures__canonical_smiles', 'similarity'])

        # Create similar drug objects and add ChEMBL info
        sim_drugs_list = [(Drug(drug["molecule_structures"]["canonical_smiles"],
                                name=drug["pref_name"],
                                chembl_id=drug["molecule_chembl_id"],
                                max_phase=drug["max_phase"],
                                indication=drug["indication_class"],
                                first_approval=drug["first_approval"]), drug["similarity"]) for drug in res]
        
        self.similarity = [(drug["molecule_chembl_id"], drug["similarity"]) for drug in res]

        return sim_drugs_list


    def smiles_similarity(mol1: str, mol2: str) -> float:
        # Create TfidfVectorizer object
        vectorizer = TfidfVectorizer()

        # Generate molecule vector
        mol1_tfidf = vectorizer.fit_transform([mol1])
        mol2_tfidf = vectorizer.fit_transform([mol2])

        # Calculate cosine similarity
        cosine_sim = cosine_similarity(mol1_tfidf, mol2_tfidf)

        return cosine_sim
    



