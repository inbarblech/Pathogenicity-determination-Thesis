from pathlib import Path
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.io import to_image
from shiny import App, Inputs, Outputs, Session, reactive, render, ui
from sklearn.decomposition import PCA
import nest_asyncio
nest_asyncio.apply()


# Function to create contour plot using PCA results and selected variants
def create_contour_plot(pca_df, selected_variants):
    # Split the PCA dataframe into Benign and Pathogenic based on 'clinsig_binary'
    Benign = pca_df[pca_df['clinsig_binary'] == 0]
    Pathogenic = pca_df[pca_df['clinsig_binary'] == 1]

    # Initialize a plotly figure
    fig = go.Figure()

    # Add Benign contour plot
    fig.add_trace(go.Histogram2dContour(
        x=Benign['PC1'], y=Benign['PC2'],
        colorscale='Greens', name='Benign',
        contours=dict(showlines=False), # Remove contour lines
        hoverinfo='skip' # Skip hover info for contours
    ))

    # Add Pathogenic contour plot
    fig.add_trace(go.Histogram2dContour(
        x=Pathogenic['PC1'], y=Pathogenic['PC2'],
        colorscale='Oranges', name='Pathogenic',
        contours=dict(showlines=False), # Remove contour lines
        hoverinfo='skip' # Skip hover info for contours
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
        title='Contour Plots for Boost Predictors', # Plot title
        xaxis_title='PC1', # X-axis label
        yaxis_title='PC2', # Y-axis label
        legend_title='Clinical Significance', # Legend title
        template='plotly_white' # Use a white background template
    )

    return fig

# Function to compute PCA and return a dataframe with the first two principal components
def compute_pca_df(features_df):
    # Initialize PCA to reduce dimensions to 2 components
    pca = PCA(n_components=2)
    # Fit the PCA on the features and transform to get the components
    pca_components = pca.fit_transform(features_df)
    # Create a DataFrame for the PCA components
    pca_df = pd.DataFrame(pca_components, columns=['PC1', 'PC2'])
    return pca_df

# Load the full data and predictions data
variants_file = Path(__file__).parent / "data.csv"
predictions = Path(__file__).parent / "predictions.csv"
features = ["stability_delta", "blosum", "hydrophobicity_delta", "volume_delta", "plddt_residue", "opra_delta", "oda_delta",
            "sasa_delta", "RSA_WT"]
full_df = pd.read_csv(variants_file)
predictions_df = pd.read_csv(predictions)

# Select specific features for PCA from the full dataset
features_df = full_df[features]

# Compute PCA for the selected features
pca_df = compute_pca_df(features_df)
# Add clinical significance and variant IDs to the PCA dataframe
pca_df['clinsig_binary'] = full_df['pathogenicity'].map({'benign': 0, 'pathogenic': 1})
pca_df['variant'] = full_df['variant']

# Define the user interface of the Shiny app
app_ui = ui.page_fluid(
    ui.layout_sidebar(
        ui.sidebar(
            ui.input_selectize("selected_variants", "Select one or more variants:", pca_df['variant'].unique().tolist(), multiple=True),
        ),
        # Instead of using ui.panel_main(), directly pass the main panel content here
        ui.output_ui("contour_plot"),
        ui.output_table("predictors_table"),
    )
)


# Define the server logic of the Shiny app
def server(input: Inputs, output: Outputs, session: Session):
    @reactive.Calc
    def selected_variants():
        # Return the rows corresponding to the selected variant IDs
        return pca_df[pca_df['variant'].isin(input.selected_variants())]

    @output
    @render.ui
    def contour_plot():
        # Create the contour plot using the selected variants
        return create_contour_plot(pca_df, selected_variants())

    @output
    @render.table
    def predictors_table():
        # Filter the predictions dataframe based on selected variant IDs
        selected_predictions = predictions_df[predictions_df['variant'].isin(input.selected_variants())]
        # Select the relevant columns and rename for better presentation
        selected_predictions = selected_predictions[[
            'variant', 'prediction'
        ]]
        selected_predictions.columns = [
            'Variant ID', 'prediction'
        ]
        return selected_predictions

# Create the Shiny app
app = App(app_ui, server)
# Run the Shiny app
app.run()