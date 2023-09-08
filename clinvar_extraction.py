import os

import pandas as pd
import general_tools as tools


def create_list_of_variant_dictionaries(txt_file):
    """Creates a list of dictionaries, each dictionary represents a variant.
    Goes over the file line by line and creates a dictionary for each variant.
    Each variant is separated by a blank line.
    Each variant dictionary contains the following keys:
    'Gene(s)', 'Protein change', 'Clinical significance', 'Review status'
    Args:
        txt_file (str): The path to the ClinVar txt file.
    Returns:
        list: A list of dictionaries, each dictionary represents a variant.
    """
    # create an empty list to hold the variants
    variants = []
    # create an empty dictionary to hold the variant
    variant = {}
    with open(txt_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            # if the line is empty, create a new dictionary
            if not line.strip():
                variants.append(variant)
                variant = {}
            # if the line is not empty, add the line to the dictionary
            else:
                # split the line to key and value
                key, value = line.split(':', 1)
                # remove the \n from the value
                value = value.strip()
                # add the key and value to the dictionary
                variant[key] = value
    # delete empty dictionaries
    variants = [variant for variant in variants if variant != {}]
    return variants


def get_dataframe_from_list_of_dict(list_of_dict) -> pd.DataFrame():
    """Creates a dataframe from a list of dictionaries.
    Args:
        list_of_dict (list): A list of dictionaries.
    returns:
        DataFrame: A dataframe created from the list of dictionaries.
    """
    # create a list of the keys in the dictionaries
    keys = list(list_of_dict[0].keys())
    # create a list of the values in the dictionaries
    values = []
    for variant in list_of_dict:
        values.append(list(variant.values()))
    # create a dataframe from the keys and values
    df = pd.DataFrame(values, columns=keys)
    # change the name of the "1. Name" column to "Name"
    df.rename(columns={'1. Name': 'Name'}, inplace=True)
    return df


def clean_variants_dataframe(df, list_of_genes):
    """Cleans the variants dataframe.
    For every row in the df: if the 'clinical significance' contains 'likely benign', delete the row.
    If the 'Review status' contains 'no assertion criteria provided' or 'conflicting interpretations of pathogenicity',
    delete the row.
    If the 'Protein change' contains more than one amino acid, split the row to two rows.
    If the 'Name' contains 'del' or 'ins', delete the row.

    Args:
        df (DataFrame): A dataframe with variants data.
    returns:
        DataFrame: A cleaned dataframe.
    """
    # for every row in the df: if the 'clinical significance' contains 'likely benign', delete the row
    df = df[~df['Clinical significance'].str.contains('Likely benign')]
    df = df[~df['Clinical significance'].str.contains('Likely pathogenic')]
    df = df[~df['Clinical significance'].str.contains('Uncertain significance')]
    # if the 'Review status' contains 'no assertion criteria provided' or 'conflicting interpretations of pathogenicity',
    # delete the row
    df = df[~df['Review status'].str.contains('no assertion criteria provided')]
    df = df[~df['Review status'].str.contains('conflicting interpretations of pathogenicity')]
    # if the 'Name' contains 'del' or 'ins', delete the row
    # if the 'Name' contains 'del' or 'ins', delete the row
    df = df[~df['Name'].str.contains('del')]
    df = df[~df['Name'].str.contains('ins')]
    # if the 'Gene' contains more than one gene, delete the row
    df = df[~df['Gene(s)'].str.contains(',')]
    # if the 'Protein change' contains more than one amino acid, split the row to two rows
    # create a new dataframe to hold the new rows
    new_df = pd.DataFrame(columns=df.columns)
    for index, row in df.iterrows():
        # if the row contains more than one variant
        if len(row['Protein change'].split(', ')) > 1:
            # split the row to two rows
            for variant in row['Protein change'].split(', '):
                # create a new row
                new_row = row.copy()
                # change the amino acid in the new row
                new_row['Protein change'] = variant
                # add the new row to the new dataframe
                new_df = new_df.append(new_row)
            # delete the row from the original dataframe
            df.drop(index, inplace=True)
    # add the new dataframe to the original dataframe
    df = df.append(new_df)

    # Creates a new "pathogenicity" column, which contains the pathogenicity of the variant
    df['Pathogenicity'] = df['Clinical significance'].apply(lambda x: 'Pathogenic' if 'Pathogenic' in x else 'Benign')

    df.drop(columns=['Accession', 'ID', 'Canonical SPDI', 'Chromosome', 'Clinical significance',
                     'Review status', 'Name', 'Location  (GRCh38)', 'Condition(s)'], inplace=True)

    # Delete rows with missing values in "Gene(s)" or "Protein change" columns
    df.dropna(subset=['Gene(s)', 'Protein change'], inplace=True)

    # Delete rows from genes that are not in list_of_genes
    df = df[df['Gene(s)'].isin(list_of_genes)]

    # Reset the index of the dataframe after dropping rows with missing values
    df.reset_index(drop=True, inplace=True)
    return df


def get_set_of_genes(path_to_gene_folders):
    """Make set of genes from the folder names in the given folder."""
    gene_names = []
    os.chdir(path_to_gene_folders)
    gene_folders = os.listdir(path_to_gene_folders)
    gene_folder_paths = [os.path.join(path_to_gene_folders, folder) for folder in gene_folders]
    for gene_folder in gene_folder_paths:
        gene_name = tools.get_folder_name_from_path_windows(gene_folder)
        if gene_name not in gene_names:
            gene_names.append(gene_name)
    return gene_names


def combine_dataframes(df1, df2) -> pd.DataFrame:
    """Combine two dataframes to one"""
    if not df1.columns.equals(df2.columns):
        raise ValueError("Columns in the two datasets are not identical.")
    # Combine the DataFrames
    combined_df = pd.concat([df1, df2], ignore_index=True)
    return combined_df


def remove_duplicate_rows(df, column1, column2):
    # Check if the specified columns exist in the DataFrame
    if column1 not in df.columns or column2 not in df.columns:
        raise ValueError("Specified columns do not exist in the DataFrame.")
    # Remove duplicate rows based on the specified columns
    df.drop_duplicates(subset=[column1, column2], keep='first', inplace=True)
    # Reset the index after removing rows
    df.reset_index(drop=True, inplace=True)
    return df


def find_identical_rows(df1, df2, column1, column2):
    # Check if the specified columns exist in both DataFrames
    if column1 not in df1.columns or column2 not in df1.columns or \
       column1 not in df2.columns or column2 not in df2.columns:
        raise ValueError("Specified columns do not exist in one or both DataFrames.")

    # Find rows with identical values in the specified columns
    merged_df = pd.merge(df1, df2, on=[column1, column2], how='inner')
    return merged_df


def remove_rows_by_values(data, remove_rows, column1, column2):
    # Check if the specified columns exist in both DataFrames
    if column1 not in data.columns or column2 not in data.columns or \
       column1 not in remove_rows.columns or column2 not in remove_rows.columns:
        raise ValueError("Specified columns do not exist in one or both DataFrames.")

    # Create a set of tuples from the values in the specified columns of "remove_rows"
    values_to_remove = set(zip(remove_rows[column1], remove_rows[column2]))

    # Filter the "data" DataFrame to remove rows with matching values
    data = data[~data.apply(lambda row: (row[column1], row[column2]) in values_to_remove, axis=1)]

    return data


if __name__ == "__main__":
    # # Create df from ClinVar txt file
    # path = 'C:\\Users\\InbarBlech\\Downloads\\clinvar_benign_missense_multiple_hearing_loss.txt'
    # paragraphs = create_list_of_variant_dictionaries(path)
    # df_benign_hearing_loss = get_dataframe_from_list_of_dict(paragraphs)
    # path = 'C:\\Users\\InbarBlech\\Downloads\\clinvar_pathogenic_missense_multiple_hearing_loss.txt'
    # paragraphs = create_list_of_variant_dictionaries(path)
    # df_pathogenic_hearing_loss = get_dataframe_from_list_of_dict(paragraphs)
    # path = 'C:\\Users\\InbarBlech\\Downloads\\clinvar_benign_missense_multiple_deafness.txt'
    # paragraphs = create_list_of_variant_dictionaries(path)
    # df_benign_deafness = get_dataframe_from_list_of_dict(paragraphs)
    # path = 'C:\\Users\\InbarBlech\\Downloads\\clinvar_pathogenic_missense_multiple_deafness.txt'
    # paragraphs = create_list_of_variant_dictionaries(path)
    # df_pathogenic_deafness = get_dataframe_from_list_of_dict(paragraphs)
    #
    # # Combine dataframes
    # df_deafness = combine_dataframes(df_pathogenic_deafness, df_benign_deafness)
    # df_hearing_loss = combine_dataframes(df_pathogenic_hearing_loss, df_benign_hearing_loss)
    # combined_df = combine_dataframes(df_hearing_loss, df_deafness)
    #
    # set_of_hearing_loss_genes_pathogenic = get_set_of_genes("C:\\Users\\InbarBlech\\OneDrive - mail.tau.ac.il\\Documents\\variants\\Pathogenic\\")
    # set_of_hearing_loss_genes_benign = get_set_of_genes("C:\\Users\\InbarBlech\\OneDrive - mail.tau.ac.il\\Documents\\variants\\Benign\\")
    # hearing_loss_genes = set_of_hearing_loss_genes_pathogenic + set_of_hearing_loss_genes_benign
    # hearing_loss_genes = set(hearing_loss_genes)
    #
    # df1 = clean_variants_dataframe(combined_df, hearing_loss_genes)
    #
    # # Remove duplicates from dataframe
    # df = remove_duplicate_rows(df1, "Gene(s)", "Protein change")
    # df.to_csv("C:\\Users\\InbarBlech\\Downloads\\clinvar_variants.csv")


    # Check duplicates in variants from dvd and clinvar

    # Create dvd file with gene and variant
    df_dvd = pd.read_csv("C:\\Users\\InbarBlech\\OneDrive - mail.tau.ac.il\\Documents\\Thesis\\Findings\\features.csv")
    df_clinvar = pd.read_csv("C:\\Users\\InbarBlech\\Downloads\\clinvar_variants.csv")
    duplicate_rows = find_identical_rows(df_dvd, df_clinvar, "gene", "variant")
    df_without_duplicates = remove_rows_by_values(df_clinvar, duplicate_rows, "gene", "variant")
    df_without_duplicates.head()
    df_without_duplicates.to_csv("C:\\Users\\InbarBlech\\Downloads\\clinvar_variants_without_duplicates_of_dvd.csv")



