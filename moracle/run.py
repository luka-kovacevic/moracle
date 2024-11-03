import streamlit as st
import pandas as pd
import numpy as np
import streamlit_ketcher as sk

from drugcomp import Drug
from prob_success import compute_prob_clin_success

# df_filename = "data/smiles_test.csv"
df_filename = "data/belka_smiles_test.csv"

# Initialize session states
if "failed_retrieval" not in st.session_state:
    st.session_state.failed_retrieval = False
if "df" not in st.session_state:
    # st.session_state.df = pd.read_csv(df_filename, sep="\t", header=0)
    st.session_state.df = pd.read_csv(df_filename, sep=",", header=0)
if "selected_chemical" not in st.session_state:
    st.session_state.selected_chemical = st.session_state.df.iloc[0]
if "selected_chemical1" not in st.session_state:
    st.session_state.selected_chemical1 = st.session_state.df.iloc[0]
if "count" not in st.session_state:
    st.session_state.count = 0

st.set_page_config(layout="wide")
st.header("MOracle 🧪🔮: Explainable clinical viability scoring in drug discovery.")
# st.write(
#     """Steps: \n
#     1. Select a candidate molecule from the DEL library.
#     2. Select a reference"""
# )

layout_mode = st.selectbox("Select Mode", options=["Single Molecule", "Comparative"])

# Helper function to compute similarity data
def get_similar_drugs_data(smiles):
    try:
        curr_drug = Drug(smiles)
        similar_drugs = curr_drug.get_similar_drugs()
        st.session_state.failed_retrieval = False
        df = pd.DataFrame.from_dict([{
            "Name": drug.name,
            "Smiles": drug.smiles,
            "Sim. Score": np.round(float(sim), 2),
            "ID": drug.chembl_id,
            "Trial Phase": drug.max_phase,
            "Ind.": drug.indication
        } for drug, sim in similar_drugs])
        return df[~df["Ind."].isna()].sort_values(by="Sim. Score", ascending=False)
    except Exception:
        st.session_state.failed_retrieval = True
        return pd.DataFrame(columns=["Name", "Smiles", "Sim. Score", "ID", "Trial Phase", "Ind."])

# Helper function to display molecule data
def display_molecule_data(label, chemical, df_similar_drugs, prob_success=None, binding_probs=None):
    chemical_name = chemical["Name"]
    molecule_smiles = chemical["Smiles"]
    data_text = f"⚛️ Current Chemical: {chemical_name}\n🗃️ SMILES: {molecule_smiles}"
    
    if binding_probs:
        data_text += "\n\nProbability of binding affinity:\n" + "\n".join([f"({prot}): {val}" for prot, val in binding_probs.items()])
    else:
        data_text += "\n\nProbability of binding affinity:\n (BRD4): n/a\n (ALB): n/a\n (FINAL): n/a"
    
    data_text += f"\n\n Historical clinical viability score : {prob_success or 'n/a'}"
    st.text_area(label=label, value=data_text, height=240)

# Helper function to add new molecule
def add_new_molecule(new_smiles):
    st.session_state.count += 1
    new_row = pd.DataFrame({"Name": f"new_mol{st.session_state.count}", "Smiles": new_smiles,}, index=[0])
    st.session_state.df = pd.concat([new_row, st.session_state.df], ignore_index=True)
    st.success(f"New molecule '{new_row.iloc[0]['Name']}' added successfully.")
    st.rerun()  # Refresh the app to display updated DataFrame

if layout_mode == "Single Molecule":

    st.subheader("Select molecule:")
    event1 = st.dataframe(st.session_state.df, height=400, key="data1", on_select="rerun", selection_mode="single-row", column_order=["Name", "Smiles"])

    if "selected_index1" not in st.session_state:
        st.session_state.selected_index1 = -1  # Set a default value if not initialized

    if len(event1.selection["rows"]) > 0:
        selected_index1 = event1.selection["rows"][0]
        if selected_index1 != st.session_state.selected_index1:
            st.session_state.selected_index1 = selected_index1
            st.session_state.selected_chemical = st.session_state.df.iloc[selected_index1]
            st.rerun() 

    chemical = st.session_state.selected_chemical
    similar_drugs_df = get_similar_drugs_data(chemical["Smiles"])
    prob_success = np.round(compute_prob_clin_success(similar_drugs_df), 2)
    binding_probs = {"BRD4": 0.45, "ALB": 0.67, "FINAL": 0.25}  # example values

    display_molecule_data("Molecule Data", chemical, similar_drugs_df, prob_success, binding_probs)

    st.session_state.main_chemical = chemical

    new_molecule_smiles = sk.st_ketcher(chemical["Smiles"], key=str(chemical["Name"]) + '_')
    if st.button("Add New Molecule from Visualiser 1"):
        add_new_molecule(new_molecule_smiles)

else:

    # UI setup for columns
    col1, col2 = st.columns([2, 2], vertical_alignment="bottom")

    #################
    ## Visualiser 1
    #################
    with col1:
        st.subheader("Select molecule:")
        event1 = st.dataframe(st.session_state.df, height=400, key="data1", on_select="rerun", selection_mode="single-row", column_order=["Name", "Smiles"])

        if "selected_index1" not in st.session_state:
            st.session_state.selected_index1 = -1  # Set a default value if not initialized

        if len(event1.selection["rows"]) > 0:
            selected_index1 = event1.selection["rows"][0]
            if selected_index1 != st.session_state.selected_index1:
                st.session_state.selected_index1 = selected_index1
                st.session_state.selected_chemical = st.session_state.df.iloc[selected_index1]
                st.rerun() 

        chemical = st.session_state.selected_chemical
        similar_drugs_df = get_similar_drugs_data(chemical["Smiles"])
        prob_success = np.round(compute_prob_clin_success(similar_drugs_df), 2)
        binding_probs = {"BRD4": 0.45, "ALB": 0.67, "FINAL": 0.25}  # example values

        display_molecule_data("Molecule Data", chemical, similar_drugs_df, prob_success, binding_probs)

        st.session_state.main_chemical = chemical

        new_molecule_smiles = sk.st_ketcher(chemical["Smiles"], key=str(chemical["Name"]) + '_')
        if st.button("Add New Molecule from Visualiser 1"):
            add_new_molecule(new_molecule_smiles)
            st.rerun()

    ################
    ## Visualiser 2
    ################
    with col2:
        st.subheader("Select reference molecule:")
        # st.session_state.option = st.selectbox("Molecule Source:", ("Belka", "ChEMBL"))
        # Display selectbox with the session state value as the default
        if "selected_source_option" not in st.session_state:
            st.session_state.selected_source_option = "Belka"  # Set a default option if not initialized

        # Now that we know `selected_source_option` exists, we can safely run the selectbox
        st.session_state.selected_source_option = st.selectbox(
            "Molecule Source:",
            options=["Belka", "ChEMBL"],
            index=["Belka", "ChEMBL"].index(st.session_state.selected_source_option)  # Set the initial option
        )
        
        if st.session_state.selected_source_option == "Belka":
            event2 = st.dataframe(st.session_state.df, height=300, key=str(st.session_state.main_chemical) + "2", on_select="rerun", selection_mode="single-row", column_order=["Name", "Smiles"])

            if len(event2.selection["rows"]) > 0:
                st.session_state.selected_chemical1 = st.session_state.df.iloc[event2.selection["rows"][0]]
            chemical = st.session_state.selected_chemical1
            similar_drugs_df = get_similar_drugs_data(chemical["Smiles"])
            prob_success = np.round(compute_prob_clin_success(similar_drugs_df), 2)
        else:
            chemical = st.session_state.selected_chemical1  
            st.session_state.similar_drugs_df = get_similar_drugs_data(st.session_state.main_chemical["Smiles"])
            prob_success = 'n/a'
            event3 = st.dataframe(st.session_state.similar_drugs_df, height=300, key=str(st.session_state.main_chemical)+ str(chemical["Name"])+"_2", on_select="rerun", selection_mode="single-row")

            if len(event3.selection["rows"]) > 0:
                st.session_state.selected_chemical1 = st.session_state.similar_drugs_df.iloc[event3.selection["rows"][0]]
            prob_success = "n/a"


        display_molecule_data("Reference Molecule Data", chemical, None, prob_success)
        new_molecule_smiles = sk.st_ketcher(chemical["Smiles"], key= str(chemical['Name'])+'_1')
        
        if st.button("Add New Molecule from Visualiser 2"):
            add_new_molecule(new_molecule_smiles)
            st.rerun()
