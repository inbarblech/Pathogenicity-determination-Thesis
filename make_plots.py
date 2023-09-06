import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt



def get_percentage_by_pathogenicity(data, feature):
    """Get the percentage of each instance of feature, in each pathogenicity group
    Example: feature = 'secondary_structure' ->
    pathogenicity  secondary_structure percentage
    benign         alpha                 0.000000
                    beta                  0.000000
                     loop                  0.000000
                    turn                  0.000000
    pathogenic     alpha                 0.000000
                    beta                  0.000000
                    loop                  0.000000
                    turn                  0.000000

    Usage: get_percentage_by_pathogenicity(data_file, 'secondary_structure')
    """
    df = pd.read_csv(data)
    # split the dataset according to pathogenicity
    df_pathogenic = df[df['pathogenicity'] == "pathogenic"]
    df_non_pathogenic = df[df['pathogenicity'] == "benign"]

    # Get the percentage of each instance of feature, in each pathogenicity group
    a = df_pathogenic.groupby(feature).count() / len(df_pathogenic)
    b = df_non_pathogenic.groupby(feature).count() / len(df_non_pathogenic)
    print(a)
    print(b)



def create_diverged_bar_plot(df):
    """ Call:
    features_df = pd.read_csv(path_to_csv_of_features, header=0)
    df = get_dataframe_of_number_of_variants_per_gene_per_pathogenicity(features_df)
    plots.create_diverged_bar_plot(df)"""
    # Sort the DataFrame by the 'gene' column to ensure the same order in both graphs
    # df_sorted = df.sort_values(by='gene', ascending=True)

    font_color = '#525252'
    hfont = {'fontname': 'David'}
    facecolor = '#eaeaf2'
    index = df.index
    column0 = df['pathogenic']
    column1 = df['benign']
    title0 = 'Pathognic variants per gene'
    title1 = 'Benign variants per gene'

    fig, axes = plt.subplots(figsize=(10, 5), facecolor=facecolor, ncols=2)
    fig.tight_layout()
    axes[0].barh(index, column0, align='center', color='red', zorder=10)
    axes[0].set_title(title0, fontsize=18, pad=15, color='red', **hfont)
    axes[1].barh(index, column1, align='center', color='lime', zorder=10)
    axes[1].set_title(title1, fontsize=18, pad=15, color='lime', **hfont)

    # If you have positive numbers and want to invert the x-axis of the left plot
    axes[0].invert_xaxis()
    axes[0].set_yticks([])  # Remove y ticks from the left plot
    axes[1].set(yticks=df.index, yticklabels=df['gene'])

    axes[1].tick_params(left=False)

    axes[1].set_xticks([10, 50, 400])
    axes[1].set_xticklabels([10, 50, 400])

    for label in axes[0].get_yticklabels():
        label.set(fontsize=5, color=font_color, **hfont)
    for label in axes[0].get_xticklabels():
        label.set(fontsize=10, color=font_color, **hfont)
    for label in axes[1].get_yticklabels():
        label.set(fontsize=2, color=font_color, **hfont)
    for label in axes[1].get_xticklabels():
        label.set(fontsize=10, color=font_color, **hfont)

    plt.subplots_adjust(wspace=0, top=0.85, bottom=0.1, left=0.18, right=0.95)
    filename = 'genes_bidirectional'
    plt.savefig(f"C:\\Users\\InbarBlech\\Downloads\\{filename}.png", dpi=2000, facecolor=facecolor)
    plt.show()


if __name__ == "__main__":
    data_file = "C:\\Users\\InbarBlech\\OneDrive - mail.tau.ac.il\\Documents\\Thesis\\Findings\\features.csv"
    # Transform to data frame
    df = pd.read_csv(data_file)
    df_pathogenic = df[df['pathogenicity'] == "pathogenic"]
    df_non_pathogenic = df[df['pathogenicity'] == "benign"]