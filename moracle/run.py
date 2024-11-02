import streamlit as st
import pandas as pd
import streamlit_ketcher as sk

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

st.session_state.df = smiles_df

col1, col2 = st.columns([2, 2])

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
            st.session_state.df,
            height=300,
            key="data1",
            on_select="rerun",
            selection_mode="single-row",
            column_order=["name", "smiles"],
        )

        st.subheader("Visualiser (Mol. 1)")

        if len(event1.selection["rows"]) > 0:
            st.session_state.selected_chemical = smiles_df.iloc[event1.selection["rows"][0]]
        else:
            st.session_state.selected_chemical = smiles_df.iloc[0]

        chemical_name = st.session_state.selected_chemical["name"]
        molecule_smiles = st.session_state.selected_chemical["smiles"]

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

        event2 = st.dataframe(
            st.session_state.df,
            height=300,
            key="data2",
            on_select="rerun",
            selection_mode="single-row",
            column_order=["name", "smiles"],
        )

        st.subheader("Visualiser (Mol. 2)")

        if len(event2.selection["rows"]) > 0:
            st.session_state.selected_chemical1 = smiles_df.iloc[event2.selection["rows"][0]]
        else: 
            st.session_state.selected_chemical = smiles_df.iloc[0]

        chemical_name1 = st.session_state.selected_chemical1["name"]
        molecule_smiles1 = st.session_state.selected_chemical1["smiles"]

        sk.st_ketcher(molecule_smiles1, key=chemical_name1 + '_1')