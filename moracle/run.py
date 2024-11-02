import streamlit as st
import pandas as pd
import streamlit_ketcher as sk

st.set_page_config(layout="wide")

smiles_df = pd.read_csv("moracle/data/smiles_test.csv", sep="\t", header=0)

col1, col2 = st.columns([4, 4])

st.header("🧪 MOracle")

st.write(
    """Molecule Oracle (MOracle) a tool for understanding differences in
         molecule-protein binding affinity between molecule structures in
         high-throughput screens."""
)

if "selected_index" not in st.session_state:
    st.session_state.selected_chemical = smiles_df.iloc[0]

st.session_state.df = smiles_df

st.subheader("Select molecule:")

st.write("Scroll through the list and select a chemical.")

event = st.dataframe(
    st.session_state.df,
    height=300,
    key="data",
    on_select="rerun",
    selection_mode=["multi-row", "multi-column"],
)

st.subheader("Visualiser")

if len(event.selection["rows"]) > 0:
    st.session_state.selected_chemical = smiles_df.iloc[event.selection["rows"][0]]

chemical_name = st.session_state.selected_chemical["name"]
molecule_smiles = st.session_state.selected_chemical["smiles"]

sk.st_ketcher(molecule_smiles)
