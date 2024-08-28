import streamlit as st
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import OneHotEncoder
import plotly.graph_objects as go


st.set_page_config(
    page_title="PredHL: Analysis",
    page_icon=":dna:",
    layout="wide",
    initial_sidebar_state="expanded")


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

    # Feature that requires one-hot encoding  #TODO: include secondary_structure in features and encode it
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
    plt.scatter(selected_point['PC1'], selected_point['PC2'], color='blue',
                label=f'{selected_variant} ({selected_gene})', )
    plt.text(selected_point['PC1'], selected_point['PC2'], f'{selected_variant}', fontsize=12,
             ha='right')

    # Plot labels and title
    plt.title('PCA of Variant Spectrum with Pathogenic/Benign Regions')
    plt.legend()

    # Remove numbers on axes
    plt.xticks([])
    plt.yticks([])

    # Show plot in Streamlit
    st.pyplot(plt)


def compute_pca_df(features_df):
    # Initialize PCA to reduce dimensions to 2 components
    pca = PCA(n_components=2)
    # Fit the PCA on the features and transform to get the components
    pca_components = pca.fit_transform(features_df)
    # Create a DataFrame for the PCA components
    pca_df = pd.DataFrame(pca_components, columns=['PC1', 'PC2'])
    return pca_df


def create_contour_plot(pca_df, selected_variants):
    # Split the PCA dataframe into Benign and Pathogenic based on 'clinsig_binary'
    Benign = pca_df[pca_df['clinsig_binary'] == 0]
    Pathogenic = pca_df[pca_df['clinsig_binary'] == 1]

    # Initialize a plotly figure
    fig = go.Figure()

    # Add Pathogenic contour plot
    fig.add_trace(go.Histogram2dContour(
        x=Pathogenic['PC1'], y=Pathogenic['PC2'],
        colorscale='Oranges', name='Pathogenic',
        contours=dict(showlines=False),  # Remove contour lines
        hoverinfo='skip'  # Skip hover info for contours
    ))

    # Add Benign contour plot
    fig.add_trace(go.Histogram2dContour(
        x=Benign['PC1'], y=Benign['PC2'],
        colorscale='Greens', name='Benign',
        contours=dict(showlines=False),  # Remove contour lines
        hoverinfo='skip'  # Skip hover info for contours
    ))

    # Add scatter plot for selected variants
    if not selected_variants.empty:
        fig.add_trace(go.Scatter(
            x=selected_variants['PC1'], y=selected_variants['PC2'],
            mode='markers+text', text=selected_variants['variant'],
            marker=dict(color='gray', size=10), # Gray color for selected variants
            textposition='top center' # Text position for variant IDs
        ))

    # Update the layout of the figure
    fig.update_layout(
        xaxis_title='PC1', # X-axis label
        yaxis_title='PC2', # Y-axis label
        legend_title='Clinical Significance', # Legend title
        template='plotly_white' # Use a white background template
    )

    return fig


# Main function to run the Streamlit app
def main():
    st.write("## Variant Spectrum Visualization")

    # Define the paths to your files
    data_file = 'data/data.csv'  # Replace with your relative path
    predictions_file = 'data/predictions.csv'  # Replace with your relative path
    features = ["stability_delta", "blosum", "hydrophobicity_delta", "volume_delta", "plddt_residue", "opra_delta",
                "oda_delta",
                "sasa_delta", "RSA_WT"]

    # Load data
    data = pd.read_csv(data_file)
    # data = load_data(data_file, predictions_file)

    # Load predictions
    predictions_df = pd.read_csv(predictions_file)


    # Select gene and variant from user input
    selected_gene = st.selectbox('Select a gene to filter:', data['gene'].unique())
    filtered_data = data[data['gene'] == selected_gene]
    selected_variant = st.multiselect('Select a variant to highlight:', options=filtered_data['variant'].unique())

    # Select specific features for PCA from the full dataset
    features_df = filtered_data[features]

    pca_df = compute_pca_df(features_df)
    # Add clinical significance and variant IDs to the PCA dataframe
    pca_df['clinsig_binary'] = filtered_data['pathogenicity'].map({'benign': 0, 'pathogenic': 1})
    pca_df['variant'] = filtered_data['variant']

    # Filter PCA dataframe based on selected variants
    if selected_variant:
        selected_df = pca_df[pca_df['variant'].isin(selected_variant)]
    else:
        selected_df = pd.DataFrame(columns=['PC1', 'PC2', 'variant'])

    # Create and display the contour plot
    fig = create_contour_plot(pca_df, selected_df)
    st.plotly_chart(fig)

    # Display the predictions table
    st.text("Predictions for Selected Variants")
    if selected_variant:
        selected_predictions = predictions_df[predictions_df['variant'].isin(selected_variant)]
        # map pathogenicity from 0 and 1 to benign and pathogenic
        selected_predictions['prediction'] = selected_predictions['prediction'].map({0: 'Benign', 1: 'Pathogenic'})
        selected_predictions = selected_predictions[['variant', 'prediction']]
        selected_predictions.columns = ['Variant', 'Prediction']
        st.dataframe(selected_predictions, hide_index=True)


# Run the app
if __name__ == "__main__":
    main()
