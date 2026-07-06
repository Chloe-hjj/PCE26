import streamlit as st
import pandas as pd

st.title("Residue Library")
#st.header("Welcome to VLP OSS")

st.subheader("Lysine")
df = pd.read_csv('data/scaffolds/lysine-choices.csv')

unique_scaffolds = df['scaffold'].unique()
option = st.selectbox("selection", unique_scaffolds)
if option:
    option_location = st.selectbox("location", ['inner', 'outer'])
    if option_location:
        st.dataframe(df.query('scaffold == @option').query('location == @option_location'))


st.subheader("The Candidate LYS")
df2 = pd.read_csv('data/scaffolds/lysine-candidates.csv')
st.dataframe(df2)

st.subheader("Aspartate/ Glutamate")
df3 = pd.read_csv('data/poi/base-choices.csv')

unique_poi = df3['poi'].unique()
option = st.selectbox("selection", unique_poi)

st.dataframe(df3.query('poi == @option'))


st.subheader("The Candidate ASP/GLU")
df4 = pd.read_csv('data/poi/best-bases.csv')
st.dataframe(df4)