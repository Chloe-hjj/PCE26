import streamlit as st
import pandas as pd

st.title("Generating Boltz2 Input")
st.write("Here you can upload .pdb / .cif files to generate boltz-2 input.")

pdb_file1 = st.file_uploader("Upload a scaffold structure", type = ["pdb", "cif"])
pdb_file2 = st.file_uploader("Upload a structure of protein of interest", type = ["pdb", "cif"])

if pdb_file1 is not None and pdb_file2 is not None:
    st.info(f"Loaded: {pdb_file1.name}")
    st.info(f"Loaded: {pdb_file2.name}")
    
    pdb1_sfring = pdb_file1.getvalue().decode("utf-8")
    pdb2_sfring = pdb_file2.getvalue().decode("utf-8")
    
    st.success("File ready fpr processing!")
    
