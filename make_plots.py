import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import general_tools as tools
from scipy import stats


def make_density_plot_per_feature_per_group(data, feature, title=None):
    """Usage:
    data_file = "C:\\Users\\InbarBlech\\OneDrive - mail.tau.ac.il\\Documents\\Thesis\\Findings\\features.csv"
    df = pd.read_csv(data_file)
    # dataframe of first group (first graph)
    df_transmembranal = df[df["protein_contain_transmembrane"] == True]
    # dataframe of second group (second graph)
    df_globular = df[df["protein_contain_transmembrane"] == False]
    make_density_plot_per_feature_per_group(df_transmembranal, "hydrophobicity_MUT", "hydrophobicity_MUT for transmembranal protein")
    make_density_plot_per_feature_per_group(df_globular, "hydrophobicity_MUT", "hydrophobicity_MUT for globular protein")"""
    df1, df2 = tools.split_dataframe_by_group(data, "pathogenicity", "benign", "pathogenic")

    # Calculate descriptive statistics for each group
    mean1, mean2 = df1[feature].mean(), df2[feature].mean()
    std1, std2 = df1[feature].std(), df2[feature].std()
    median1, median2 = df1[feature].median(), df2[feature].median()

    # Create density plots for each sub-DataFrame
    plt.figure(figsize=(12, 6))
    sns.kdeplot(data=df1[feature], label=f'{"Benign"} - {feature}', shade=True, color="green")
    sns.kdeplot(data=df2[feature], label=f'{"Pathogenic"} - {feature}', shade=True, color="red")

    # Customize plot labels and title
    plt.xlabel(feature)
    plt.ylabel("Density")
    if title is None:
        plt.title(f"Density Plot for {feature}")
    else:
        plt.title(title)

    # Show legend
    plt.legend()

    # plt.savefig(f"C:\\Users\\InbarBlech\\OneDrive - mail.tau.ac.il\\Documents\\Thesis\\Findings\\Plots\\{title}.png", dpi=1200)

    # Perform t-test and calculate p-value
    t_stat, p_value = stats.ttest_ind(df1[feature].dropna(), df2[feature].dropna())
    # Format p-value with two decimal places, with a minimum of two digits after the dot
    formatted_p_value = '{:.2f}'.format(p_value) if p_value >= 0.01 else '0.00'

    # # Optional: Add the p-value to the plot
    # plt.figure(figsize=(12, 6))
    # sns.kdeplot(data=df1[feature], label=f'Benign - {feature}', shade=True, color="green")
    # sns.kdeplot(data=df2[feature], label=f'Pathogenic - {feature}', shade=True, color="red")
    # plt.xlabel(feature)
    # plt.ylabel("Density")
    # if title is None:
    #     plt.title(f"Density Plot for {feature}")
    # else:
    #     plt.title(title)
    # plt.legend()
    # plt.figtext(0.15, 0.85, f'T-test p-value: {formatted_p_value}', bbox=dict(facecolor='white', alpha=0.5))
    # plt.show()

    # Display statistics on the plot
    stats_text = (
        f"Benign: mean={mean1:.2f}, std={std1:.2f}, median={median1:.2f}\n"
        f"Pathogenic: mean={mean2:.2f}, std={std2:.2f}, median={median2:.2f}\n"
        f"T-test p-value: {formatted_p_value}"
    )
    plt.figtext(0.15, 0.75, stats_text, bbox=dict(facecolor='white', alpha=0.5))

    # Display the plot
    plt.show()


def make_unsmoothed_density_plot_per_feature_per_group(data, feature, title=None):
    """
    Create density plots for each sub-DataFrame, without smoothing.
    """

    # Create a bar plot of value frequencies
    plt.figure(figsize=(12, 6))
    # Group the data by 'feature' and 'pathogenicity' and count the occurrences
    grouped = data.groupby(feature)['pathogenicity'].value_counts().unstack().fillna(0)

    # Set the width of each bar
    bar_width = 0.35

    # Get the unique values from Column B (x-axis values)
    x_values = grouped.index

    # Get the number of unique values in Column A
    num_unique_values_a = len(data['pathogenicity'].unique())

    # Calculate the position of bars for two groups
    group1_positions = range(len(x_values))
    group2_positions = [pos + bar_width for pos in group1_positions]

    # Create the bar plot
    plt.bar(group1_positions, grouped.iloc[:, 0], width=bar_width, label='Group 1')
    plt.bar(group2_positions, grouped.iloc[:, 1], width=bar_width, label='Group 2')

    # Customize the plot
    plt.xlabel('Values in pathogenicity')
    plt.ylabel('Frequency')
    plt.title('Frequency of Values in Column B by Group')
    plt.xticks([pos + bar_width / 2 for pos in group1_positions], x_values)
    plt.legend()

    # Show the plot
    plt.tight_layout()
    plt.show()
    #
    # if title is None:
    #     plt.title(f"Grouped Bar Plot for {feature}")
    # else:
    #     plt.title(title)
    #
    # # Show legend
    # plt.legend(title="Pathogenicity")
    #
    # # Save the plot
    # plt.savefig(f"C:\\Users\\InbarBlech\\OneDrive - mail.tau.ac.il\\Documents\\Thesis\\Findings\\Plots\\{title}.png",
    #             dpi=1200)
    #
    # # Display the plot
    # plt.show()


def get_percentage_by_pathogenicity(df, feature) -> (pd.DataFrame, pd.DataFrame):
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
    # split the dataset according to pathogenicity
    df_pathogenic = df[df['pathogenicity'] == "pathogenic"]
    df_non_pathogenic = df[df['pathogenicity'] == "benign"]
    #
    # df_secondary_structure_pathogenic = df_pathogenic[df_pathogenic[feature]]
    # df_secondary_structure_benign = df_non_pathogenic[df_non_pathogenic[feature]]

    # Get the percentage of each instance of feature, in each pathogenicity group
    a = df_pathogenic.groupby(feature).size() / len(df_pathogenic)
    b = df_non_pathogenic.groupby(feature).size() / len(df_non_pathogenic)

    return a, b


def plot_two_pie_charts(series1, series2, title1, title2, path_to_save_fig) -> None:
    """This function saves two pie charts plots, derived from two pandas series.
    Args:
        series1
        series2:
        title1 (str): Title for the first plot.
        title2 (str): Title for the second plot.
        path_to_save_fig (str): Path to save figure to.
    Call:
        # Organize:
        data_file = "C:\\Users\\InbarBlech\\OneDrive - mail.tau.ac.il\\Documents\\Thesis\\Findings\\features.csv"
        df = pd.read_csv(data_file)
        # Get percentages using get_percentage_by_pathogenicity
        series1, series2 = plots.get_percentage_by_pathogenicity(df, column_to_plot_by)  # For example, secondary structure
        # Call this function:
        plots.plot_two_pie_charts(series1, series2, 'title1', 'title2', path_to_save_fig)
    """
    color_mapping = {"Loop": 'mediumpurple', "Helix": 'moccasin', "Beta strand": 'powderblue', "Turn": 'lavender'}

    labels1 = series1.index
    percentages1 = series1.iloc[:]
    plt.figure(figsize=(10, 5))  # Optional: Adjust the figure size
    plt.subplot(121)  # Create the first subplot for 'pathogenic'
    colors1 = [color_mapping[label] for label in labels1]
    plt.pie(percentages1, labels=labels1, autopct='%1.1f%%', startangle=140, colors=colors1)
    plt.title(title1, fontsize=12)
    plt.axis('equal')

    labels2 = series2.index
    percentages2 = series2.iloc[:]
    plt.subplot(122)  # Create the second subplot for 'benign'
    colors2 = [color_mapping[label] for label in labels2]
    plt.pie(percentages2, labels=labels2, autopct='%1.1f%%', startangle=140, colors=colors2)
    plt.title(title2, fontsize=12)
    plt.axis('equal')

    plt.tight_layout()
    plt.savefig(path_to_save_fig)


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


def create_distribution_of_sequence_length(df: pd.DataFrame, path_to_save):
    """Creates distribution of sequence length in dataset, diveded by number of variants per gene.
    If there are multiple variants per gene, the gene is not counted more than once.
    Assumes
    sequence_length is a column in the df.
    Call:
    features_df = pd.read_csv(path_to_csv_of_features, header=0)
    plots.create_distribution_of_sequence_length(features_df)"""
    df = df.drop_duplicates(subset=['gene'], keep='first')
    sequence_length = df['sequence_length']
    plt.figure(figsize=(10, 5))
    plt.hist(sequence_length, bins=100, color='lightblue')
    plt.title('Distribution of protein sequence length', fontsize=18)
    plt.xlabel('Sequence length', fontsize=14)
    plt.ylabel('Number of genes in dataset', fontsize=14)
    plt.savefig(f"{path_to_save}\\sequence_length_distribution.png", dpi=2000)
    plt.show()


def get_number_of_transmembrane_globular_residues_in_df(df: pd.DataFrame):
    """Returns the number of transmembrane and globular residues in the df.
    Assumes
    'transmembrane_residues' and 'globular_residues' are columns in the df.
    Call:
    features_df = pd.read_csv(path_to_csv_of_features, header=0)
    plots.get_number_of_transmembrane_globular_residues_in_df(features_df)"""
    transmembrane_residues = len(df[df['is_residue_transmembranal']])
    globular_residues = len(df[~df['is_residue_transmembranal']])
    return transmembrane_residues, globular_residues


if __name__ == "__main__":
    df = pd.read_csv('/home/inbar/results/combined_with_source.csv', header=0)
    print(len(df))
    # df = pd.read_csv(f'/home/inbar/results/gene_specific_df/{gene}.csv', header=0)
    # print(f"transmename and globular residues in gene {gene}: {get_number_of_transmembrane_globular_residues_in_df(df)}")
    #
    # Use make_density_plot_per_feature_per_group
    make_density_plot_per_feature_per_group(df, "plddt_residue", "plddt_residue for DVD variants")