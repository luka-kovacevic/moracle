import streamlit as st
import pandas as pd
import numpy as np
import streamlit_ketcher as sk
from streamlit.components.v1 import html
from mol_viewer import gen_3dmol_vis, run_wrapper

from drugcomp import Drug
from prob_success import compute_prob_clin_success

count = 0
# Function to add a new row to the DataFrame
def add_row(name, smiles, df):
    new_row = pd.DataFrame({"Name": [name], "SMILES": [smiles]})
    st.session_state.df = pd.concat([st.session_state.df, new_row], ignore_index=True)

# Initialize `failed_retrieval` in session state
if "failed_retrieval" not in st.session_state:
    st.session_state.failed_retrieval = False

st.set_page_config(layout="wide")

st.session_state.df = pd.read_csv("data/smiles_test.csv", sep="\t", header=0)

st.header("🧪 MOracle")

st.write(
    """Molecule Oracle (MOracle) a tool for understanding differences in
         molecule-protein binding affinity between molecule structures in
         high-throughput screens."""
)

if "selected_chemical" not in st.session_state:
    st.session_state.selected_chemical = st.session_state.df.iloc[0]
    st.session_state.selected_chemical1 = st.session_state.df.iloc[0]

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
            st.session_state.df,
            height=400,
            key="data1",
            on_select="rerun",
            selection_mode="single-row",
            column_order=["Name", "Smiles"],
        )

        st.subheader("Visualiser")

        if len(event1.selection["rows"]) > 0:
            st.session_state.selected_chemical = st.session_state.df.iloc[event1.selection["rows"][0]]
        else:
            st.session_state.selected_chemical = st.session_state.df.iloc[0]

        chemical_name = st.session_state.selected_chemical["Name"]
        molecule_smiles = st.session_state.selected_chemical["Smiles"]

        curr_drug = Drug(molecule_smiles)
        similar_drugs = curr_drug.get_similar_drugs()
        similar_drugs_df = pd.DataFrame.from_dict([{"Name":drug.name, "Smiles":drug.smiles, "Sim. Score": np.round(float(sim), 2), "ID":drug.chembl_id, "Trial Phase":drug.max_phase, "Ind.":drug.indication} for drug, sim in similar_drugs])
        similar_drugs_df = similar_drugs_df[~similar_drugs_df["Ind."].isna()].sort_values(by="Sim. Score", ascending=False)


        txt = st.text_area(label="Molecule Data",
            value="⚛️ Current Chemical: " + chemical_name +
            "\n🗃️ Smiles: " + molecule_smiles +
            "\n" + 
            "\nProbability of binding affinity:" + 
            "\n (BRD4): 0.45 (TEST_VAL)" +
            "\n (ALB): 0.67 (TEST_VAL)" + 
            "\n (FINAL): 0.25 (TEST_VAL)" +
            "\n" +
            "\nProbability of clinical success: " + str(np.round(compute_prob_clin_success(similar_drugs_df), 2)),
            height=240
        )

        new_molecule_smiles = sk.st_ketcher(molecule_smiles, key=chemical_name + '_')

        if new_molecule_smiles not in st.session_state.df["Smiles"]:  
                new_row = pd.DataFrame({"Name":"new_mol"+str(count), "Smiles": new_molecule_smiles, "ID": "", "Trial Phase": "", "Ind.": ""}, index=[0])
                st.session_state.df = pd.concat([st.session_state.df, new_row], ignore_index=True)
                count += 1

        html(run_wrapper("./data/complex_example.zip"), height=1000)


################
## Visualiser 2
################
with col2:
    with st.container():

        ########################
        ## Molecule selection 1
        ########################

        st.subheader("Select reference molecule:")

        st.write("Scroll through the list and select a chemical.")

        option = st.selectbox("Molecule Source:",
                             ("Belka", "ChEMBL"))


        # st.subheader("Select molecule:")

        # st.write("Scroll through the list and select a chemical.")

        if option == "Belka":

            event2 = st.dataframe(
                st.session_state.df,
                height=300,
                key="data2",
                on_select="rerun",
                selection_mode="single-row",
                column_order=["Name", "Smiles"],
            )

            if len(event2.selection["rows"]) > 0:
                st.session_state.selected_chemical1 = st.session_state.df.iloc[event2.selection["rows"][0]]
            else: 
                st.session_state.selected_chemical = st.session_state.df.iloc[0]

            chemical_name1 = st.session_state.selected_chemical1["Name"]
            molecule_smiles1 = st.session_state.selected_chemical1["Smiles"]

            st.subheader("Visualiser (Reference)")

            html(run_wrapper("./data/complex_example.zip"), height=1000)

            curr_drug = Drug(st.session_state.selected_chemical["Smiles"])

            curr_drug = Drug(molecule_smiles1)

            similar_drugs = curr_drug.get_similar_drugs()
            similar_drugs_df = pd.DataFrame.from_dict([{"Name":drug.name, "Smiles":drug.smiles, "Sim. Score": np.round(float(sim), 2), "ID":drug.chembl_id, "Trial Phase":drug.max_phase, "Ind.":drug.indication} for drug, sim in similar_drugs])
            similar_drugs_df = similar_drugs_df[~similar_drugs_df["Ind."].isna()].sort_values(by="Sim. Score", ascending=False)

            txt = st.text_area(label="Molecule Data",
                value="⚛️ Current Chemical: " + chemical_name +
                "\n🗃️ Smiles: " + molecule_smiles +
                "\n" + 
                "\nProbability of binding affinity:" + 
                "\n (BRD4): 0.45 (TEST_VAL)" +
                "\n (ALB): 0.67 (TEST_VAL)" + 
                "\n (FINAL): 0.25 (TEST_VAL)" +
                "\n" +
                "\nProbability of clinical success: " + str(np.round(compute_prob_clin_success(similar_drugs_df), 2)),
                height=240, 
                key= chemical_name + "_textelem"
            )

            # txt = st.text_area(
            #     f"⚛️ Current Chemical: {chemical_name}",
            #     (
            #         f"🗃️ Smiles: {molecule_smiles}\n"
            #         "Probability of binding affinity (BRD4): 0.45 (TEST_VAL)\n"
            #         "                                (ALB): 0.67 (TEST_VAL)\n"
            #         "                                (FINAL): 0.25 (TEST_VAL)\n"
            #         "Probability of clinical success ()"
            #     )
            # )

            new_molecule_smiles = sk.st_ketcher(molecule_smiles1, key=chemical_name1 + '_1')

            if new_molecule_smiles not in st.session_state.df["Smiles"]:
                new_row = pd.DataFrame({"Name":"new_mol"+str(count), "Smiles": new_molecule_smiles, "ID": "", "Trial Phase": "", "Ind.": ""}, index=[0])
                smiles_df = pd.concat([st.session_state.df, new_row], ignore_index=True)
                count += 1

        else: 

            try:
                # Attempt to create Drug instance and retrieve similar drugs
                curr_drug = Drug(st.session_state.selected_chemical["Smiles"])
                similar_drugs = curr_drug.get_similar_drugs()
                
                # Build DataFrame from similar drugs data
                similar_drugs_df = pd.DataFrame.from_dict([
                    {
                        "Name": drug.name,
                        "Smiles": drug.smiles,
                        "Sim. Score": np.round(float(sim), 2),
                        "ID": drug.chembl_id,
                        "Trial Phase": drug.max_phase,
                        "Ind.": drug.indication
                    }
                    for drug, sim in similar_drugs
                ])

                # Filter and sort DataFrame
                st.session_state.df_drugs = similar_drugs_df[~similar_drugs_df["Ind."].isna()].sort_values(by="Sim. Score", ascending=False)
                
                # Display DataFrame
                event2 = st.dataframe(
                    st.session_state.df_drugs,
                    height=300,
                    key="data2",
                    on_select="rerun",
                    selection_mode="single-row",
                )
                
                # Set `failed_retrieval` to False if successful
                st.session_state.failed_retrieval = False

            except Exception:
                # Handle any errors and set `failed_retrieval` to True
                st.session_state.failed_retrieval = True

            if len(event2.selection["rows"]) > 0:
                st.session_state.selected_chemical1 = st.session_state.df_drugs.iloc[event2.selection["rows"][0]]
            else: 
                st.session_state.selected_chemical1 = st.session_state.df_drugs.iloc[0]

            chemical_name1 = st.session_state.selected_chemical1["Name"]
            molecule_smiles1 = st.session_state.selected_chemical1["Smiles"]
             
            st.subheader("Visualiser (Reference)")

            st.session_state.text_content = ""

            if st.session_state.failed_retrieval:
                st.session_state.text_content = "(molecule not found!) \n"
            

            txt = st.text_area(label="Molecule Data",
                value=st.session_state.text_content + "Current Chemical: " + chemical_name1 +
                "\n SMILES: " + molecule_smiles1 +
                "\n" + 
                "\nProbability of binding affinity:" + 
                "\n (BRD4): n/a" +
                "\n (ALB): n/a" + 
                "\n (FINAL): n/a" +
                "\n" +
                "\nProbability of clinical success: n/a",

                height=240, 
                key= chemical_name1 + "_textelem"
            )

            new_molecule_smiles = sk.st_ketcher(molecule_smiles1, key=chemical_name1 + '_1')

            if new_molecule_smiles not in st.session_state.df_drugs["Smiles"]:
                new_row = pd.DataFrame({"Name":"new_mol"+str(count), "Smiles": new_molecule_smiles, "ID": "", "Trial Phase": "", "Ind.": ""}, index=[0])
                st.session_state.df_drugs = pd.concat([st.session_state.df_drugs, new_row], ignore_index=True)
                count += 1
                st.rerun()

            sk.st_ketcher(molecule_smiles1, key=chemical_name1 + '_1')

            html(run_wrapper("./data/complex_example.zip"))
