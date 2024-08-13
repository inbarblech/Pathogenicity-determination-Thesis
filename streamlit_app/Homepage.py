import streamlit as st
import numpy as np
import pandas as pd
import joblib


# Load the model
# model = joblib.load("model.joblib")

st.set_page_config(
    page_title="ML App",
    page_icon=":dna:",
    layout="wide",
    initial_sidebar_state="expanded")


# Header
with st.container():
    st.title("HearIPred - Hearing Impairment predictor and analysis tool")
    st.write("This app is a demo of a machine learning application for variant analysis.")
    st.write("[GitHub Repository](https://github.com/inbarblech/Pathogenicity-determination-Thesis)")
