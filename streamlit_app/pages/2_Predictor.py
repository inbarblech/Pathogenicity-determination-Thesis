import streamlit as st
import numpy as np
import pandas as pd
import joblib
from utils import check_variant_valid

# Load the model
# model = joblib.load("model.joblib")

genes_list = ["SLC26A4", "GJB2", "COL2A1", "COL4A5", "COL4A3", "MYO7A", "WFS1", "FGFR1"]

with st.container():
    st.write("## Predictor")
    left_column, right_column = st.columns(2)
    with left_column:
        st.subheader("Enter a variant")
        variant = st.text_input("Variant", placeholder="e.g. R75Q")
    with right_column:
        st.subheader("Select a gene")
        gene = st.selectbox("Select a gene", genes_list)
    # Add a button
    if st.button("Predict"):
        while variant == "" or gene == "":
            st.write("Please enter a variant and a gene.")
            break
        # Check that variant is in the correct format
        if not check_variant_valid(variant, gene):  # TODO: Implement this function
            st.write("Error. Please enter a valid variant.")

        if variant and gene:
            # Perform feature extraction
            st.write("Computing features for the variant...")
            # features = extract_features(variant, gene)
            # # Make predictions
            # st.write(f"Predicting the pathogenicity of the variant {variant} in gene {gene}...")
            # pred, prob = model.predict(features)
            # # Display the prediction
            # st.subheader("Prediction Results")
            # st.write(f"Prediction: {pred}")
            # st.write(f"Probability: {prob}")
            # st.progress(prob)
        else:
            st.write("Please enter valid variant and gene.")
    st.write("---")
