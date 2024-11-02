import streamlit as st
import pandas as pd
import numpy as np
import streamlit_ketcher as sk

from drugcomp import Drug

st.set_page_config(layout="wide")

smiles_df = pd.read_csv("data/smiles_test.csv", sep="\t", header=0)

st.header("🧪 MOracle")

st.write(
    """Molecule Oracle (MOracle) a tool for understanding differences in
         molecule-protein binding affinity between molecule structures in
         high-throughput screens."""
)

if "selected_chemical" not in st.session_state:
    st.session_state.selected_chemical = smiles_df.iloc[0]
    st.session_state.selected_chemical1 = smiles_df.iloc[0]

# st.session_state.df = smiles_df

col1, col2 = st.columns([2, 2], vertical_alignment="bottom")

#################
## Visualiser 1
################


with col1: 
    with st.container():

        ########################
        ## Molecule selection 1
        ########################

        st.subheader("Select molecule:")

        st.write("Scroll through the list and select a chemical.")

        event1 = st.dataframe(
            smiles_df,
            height=400,
            key="data1",
            on_select="rerun",
            selection_mode="single-row",
            column_order=["Name", "Smiles"],
        )

        st.subheader("Visualiser (Mol. 1)")

        if len(event1.selection["rows"]) > 0:
            st.session_state.selected_chemical = smiles_df.iloc[event1.selection["rows"][0]]
        else:
            st.session_state.selected_chemical = smiles_df.iloc[0]

        chemical_name = st.session_state.selected_chemical["Name"]
        molecule_smiles = st.session_state.selected_chemical["Smiles"]

        sk.st_ketcher(molecule_smiles, key=chemical_name + '_')


################
## Visualiser 2
################
with col2:
    with st.container():

        ########################
        ## Molecule selection 1
        ########################

        st.subheader("Select molecule:")

        st.write("Scroll through the list and select a chemical.")

        option = st.selectbox("Molecule Source:",
                             ("Belka", "ChEMBL"))


        # st.subheader("Select molecule:")

        # st.write("Scroll through the list and select a chemical.")

        if option == "Belka":

            event2 = st.dataframe(
                smiles_df,
                height=300,
                key="data2",
                on_select="rerun",
                selection_mode="single-row",
                column_order=["Name", "Smiles"],
            )

            if len(event2.selection["rows"]) > 0:
                st.session_state.selected_chemical1 = smiles_df.iloc[event2.selection["rows"][0]]
            else: 
                st.session_state.selected_chemical = smiles_df.iloc[0]

            chemical_name1 = st.session_state.selected_chemical1["Name"]
            molecule_smiles1 = st.session_state.selected_chemical1["Smiles"]

            st.subheader("Visualiser (Mol. 2)")
            sk.st_ketcher(molecule_smiles1, key=chemical_name1 + '_1')

        else: 
            curr_drug = Drug(st.session_state.selected_chemical["Smiles"])
            similar_drugs = curr_drug.get_similar_drugs()
            similar_drugs_df = pd.DataFrame.from_dict([{"Name":drug.name, "Smiles":drug.smiles, "Sim. Score": np.round(float(sim), 2), "ID":drug.chembl_id, "Trial Phase":drug.max_phase, "Ind.":drug.indication} for drug, sim in similar_drugs])

            similar_drugs_df = similar_drugs_df[~similar_drugs_df["Ind."].isna()].sort_values(by="Sim. Score", ascending=False)

            event2 = st.dataframe(
                similar_drugs_df,
                height=300,
                key="data2",
                on_select="rerun",
                selection_mode="single-row",
            )

            if len(event2.selection["rows"]) > 0:
                st.session_state.selected_chemical1 = similar_drugs_df.iloc[event2.selection["rows"][0]]
            else: 
                st.session_state.selected_chemical1 = similar_drugs_df.iloc[0]

            chemical_name1 = st.session_state.selected_chemical1["Name"]
            molecule_smiles1 = st.session_state.selected_chemical1["Smiles"]
             
            st.subheader("Visualiser (Mol. 2)")

            sk.st_ketcher(molecule_smiles1, key=chemical_name1 + '_1')