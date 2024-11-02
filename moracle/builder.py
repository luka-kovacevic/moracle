import streamlit as st
import pandas as pd
import streamlit_ketcher as sk

class Page:

    def __init__(self, df: pd.DataFrame, selection_default: int):
        self.dataframe = st.dataframe(
            df,
            height=300,
            key="data1",
            on_select="rerun",
            selection_mode="single-row",
            column_order=["name", "smiles"],
        )

        self.selection = selection_default
        self.chemical_name = df.iloc[self.selection]["names"]
        self.molecule_smiles = df.iloc[self.selection]["smiles"]

        self.molecule_drawer = sk.st_ketcher(self.molecule_smiles)