import streamlit as st
import sys
import os
import re

#Automate the path changing to the directory where the .csv is
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)

import subject_checker

st.set_page_config(layout="wide")
st.header("Welcome to the iGEM team search program !", divider='rainbow')

keywords = st.text_input("## Enter the keywords").split(',')

st.markdown("This tool was made by the team iGEM IONIS 2024")

if st.button("Start Research"):
    results = subject_checker.subject_finder(keywords)

    st.markdown(f"{len(results)} team found")

    if len(results) < 30 :
        for team in results : 
            team = team.strip('\n').split('\t')
            col1, col2, col3 = st.columns([2, 1.2, 10])
            col1.markdown(team[0])
            col2.markdown(team[1])
            col3.markdown(team[2])
    
    else:
        st.markdown("**Too many results to display abstract, only the team link and village will be shown**")
        for team in results : 
                team = team.strip('\n').split('\t')
                col1, col2 = st.columns([1, 1])
                col1.markdown(team[0])
                col2.markdown(team[1])
        




