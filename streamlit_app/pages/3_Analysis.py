import streamlit as st
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import OneHotEncoder


# Function to load data
def load_data(data_file, predictions_file):
    """Load the dataset and predictions from CSV files."""
    data = pd.read_csv(data_file)
    predictions = pd.read_csv(predictions_file)
    merged_data = pd.merge(data, predictions, on='variant')
    # Ensure that the 'gene' column is retained
    if 'gene' not in merged_data.columns:
        merged_data['gene'] = data['gene']
    return merged_data


# Function to perform PCA
def perform_pca(data, n_components=2):
    """Perform PCA on the dataset and return a DataFrame with principal components."""
    # Non-numeric features to exclude from PCA
    non_features = [
        'gene', 'variant', 'pathogenicity', 'protein_contain_transmembrane',
        'is_residue_transmembranal', 'source', 'aa_WT', 'aa_MUT',
        'Absolute value of delta', 'stability_WT', 'stability_MUT',
        'hydrophobicity_WT', 'hydrophobicity_MUT', 'volume_WT',
        'volume_MUT', 'opra_WT', 'opra_MUT', 'opra_delta', 'oda_WT',
        'oda_MUT', 'sasa_WT', 'sasa_MUT', 'sequence_length', 'RSA_MUT'
    ]

    features = [
        'blosum', 'plddt_residue', 'stability_delta', 'hydrophobicity_delta',
        'volume_delta', 'RSA_WT', 'oda_delta', 'sasa_delta'
    ]

    # Feature that requires one-hot encoding
    features_to_encode = ['secondary_structure']

    # # Perform one-hot encoding on categorical features
    # encoder = OneHotEncoder(sparse_output=False)
    # encoded_features = pd.DataFrame(encoder.fit_transform(data[features_to_encode]))
    # encoded_features.columns = encoder.get_feature_names_out(features_to_encode)
    #
    # # Combine the numeric and encoded features
    # data = data[features + features_to_encode]

    data_features = data[features]

    # Apply PCA on the numeric and encoded data
    pca = PCA(n_components=n_components)
    principal_components = pca.fit_transform(data_features)

    # Create a DataFrame with the PCA results
    pca_df = pd.DataFrame(data=principal_components, columns=[f'PC{i + 1}' for i in range(n_components)])
    pca_df['variant'] = data['variant']
    pca_df['gene'] = data['gene']
    pca_df['prediction'] = data['prediction']
    return pca_df


# Function to plot PCA results
def plot_pca(pca_df, selected_variant, selected_gene):
    """Plot PCA results with KDE plots, highlighting the selected variant and gene."""
    plt.figure(figsize=(10, 8))

    # KDE plot for pathogenic variants
    sns.kdeplot(x=pca_df[pca_df['prediction'] == 'Pathogenic']['PC1'],
                y=pca_df[pca_df['prediction'] == 'Pathogenic']['PC2'],
                fill=True, cmap="Oranges", label='Pathogenic')

    # KDE plot for benign variants
    sns.kdeplot(x=pca_df[pca_df['prediction'] == 'Benign']['PC1'],
                y=pca_df[pca_df['prediction'] == 'Benign']['PC2'],
                fill=True, cmap="Greens", label='Benign')

    # Highlight the selected variant and gene
    selected_point = pca_df[(pca_df['variant'] == selected_variant) & (pca_df['gene'] == selected_gene)]
    plt.scatter(selected_point['PC1'], selected_point['PC2'], color='black',
                label=f'{selected_variant} in {selected_gene}')
    plt.text(selected_point['PC1'], selected_point['PC2'], f'{selected_variant} ({selected_gene})', fontsize=12,
             ha='right')

    # Plot labels and title
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('PCA of Variant Spectrum with Pathogenic/Benign Regions')
    plt.legend()

    # Show plot in Streamlit
    st.pyplot(plt)


# Main function to run the Streamlit app
def main():
    st.write("## Variant and Gene Spectrum Visualization using PCA with Pathogenic/Benign Regions")

    # Define the paths to your files
    data_file = 'data/data.csv'  # Replace with your relative path
    predictions_file = 'data/predictions.csv'  # Replace with your relative path

    # Load data
    data = load_data(data_file, predictions_file)

    # Select gene and variant from user input
    selected_gene = st.selectbox('Select a gene to filter:', data['gene'].unique())
    filtered_data = data[data['gene'] == selected_gene]
    selected_variant = st.selectbox('Select a variant to highlight:', filtered_data['variant'].unique())

    # Perform PCA
    pca_df = perform_pca(data)

    # Plot results
    plot_pca(pca_df, selected_variant, selected_gene)


# Run the app
if __name__ == "__main__":
    main()
