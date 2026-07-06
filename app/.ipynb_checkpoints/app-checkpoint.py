import streamlit as st
import pandas as pd

st.title("VLPOSS")
st.header("Welcome to VLP OSS")

df = pd.read_csv('../data/scaffolds/lysine-choices.csv')

unique_scaffolds = df['scaffold'].unique()
option = st.selectbox("selection", unique_scaffolds)
if option:
    option_location = st.selectbox("location", ['inner', 'outter'])
    if option_location:
        st.dataframe(df.query('scaffold == @option').query('location == @option_location'))