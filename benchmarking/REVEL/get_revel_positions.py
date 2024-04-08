import os

import pandas as pd
import ensemble_api
import requests as req


def add_aa_position_to_df(df: pd.DataFrame) -> pd.DataFrame:
    """This function adds the amino acid position to the dataframe, under new column 'aa_pos'
    Method:
    Add a new column 'aa_pos' to the dataframe
    Iterate through the "pos_hg19" column, and for each increasment of 3, add the value to the new column
    """
    # Initialize variables to keep track of the last value and the incremental number
    last_value = None
    incremental_number = 0
    number_column = []

    # Iterate through the pos_hg19 column and calculate the numbers column
    for value in df['pos_hg19']:
        if last_value is None:
            last_value = value
            incremental_number = 1
        elif value - last_value >= 3:
            incremental_number += 1
            last_value = value
        else:
            pass
        number_column.append(incremental_number)

    # Add the numbers column to the DataFrame
    df['aa_pos'] = number_column
    return df


def add_aa_position_to_df_reverse(df: pd.DataFrame) -> pd.DataFrame:
    """This function adds the amino acid position to the dataframe, under new column 'aa_pos'
    Method:
    Add a new column 'aa_pos' to the dataframe
    Iterate through the "pos_hg19" column, and for each increasment of 3, add the value to the new column
    """
    # Initialize variables to keep track of the last value and the incremental number
    last_value = None
    incremental_number = 0
    number_column = []

    # Iterate through the pos_hg19 column and calculate the numbers column
    for value in df['pos_hg19']:
        if last_value is None:
            last_value = value
            incremental_number = 1
        elif last_value - value >= 3:
            incremental_number += 1
            last_value = value
        else:
            pass
        number_column.append(incremental_number)

    # Add the numbers column to the DataFrame
    df['aa_pos'] = number_column
    return df


def encode_amino_acids(df):
    df = df.sort_values(by='aa_pos').reset_index(drop=True)  # Sort the DataFrame by 'aa_pos' and reset the index
    aa_seq = ''

    for index, row in df.iterrows():
        aa_wt = row['aa_wt']
        aa_seq += aa_wt

    return aa_seq


def get_uniprot_url(gene_name) -> str:
    """Returns the URL for the Uniprot page for the given gene name."""
    url = f"https://rest.uniprot.org/uniprotkb/search?query=(gene:{gene_name})%20AND%20(taxonomy_id:9606)%20AND%20(reviewed:true)"
    return url


def get_uniprot_json(gene_name) -> dict:
    # The base URL for UniProt's search API
    url = get_uniprot_url(gene_name)
    # Make a request to the search API
    response = req.get(url)
    # Extract the JSON data from the response
    data = response.json()
    return data


def get_uniprot_id(gene_name) -> str:
    """Returns the Uniprot ID for the given gene name."""
    data = get_uniprot_json(gene_name)
    primary_accession = data['results'][0]['primaryAccession']
    return primary_accession


def get_sequence(gene_name) -> str:
    """Returns the sequence for the given Uniprot"""
    data = get_uniprot_json(gene_name)
    sequence = data['results'][0]['sequence']['value']
    return sequence


def create_genes_location_dict(genes_list):
    """This function creates a dictionary of genes and their genomic location.
    The keys are the gene names and the values are tuples of (chromosome, start, end, strand(1 or -1)"""
    genes_location_dict = {}
    for gene in genes_list:
        genes_location_dict[gene] = ensemble_api.get_genomic_location(gene)
        print(f"Added {gene} to the genes_location_dict. Genomic location: {genes_location_dict[gene]}")
    return genes_location_dict


def cut_dataframe_by_start_end_locations(df: pd.DataFrame, start: int, end: int, assembly_name) -> pd.DataFrame:
    """This function cuts the dataframe by the start and end locations"""
    if assembly_name == 'GRCh38':
        return df[(int(df['grch38_pos']) >= start) & (int(df['grch38_pos']) <= end)]
    elif assembly_name == 'GRCh37':
        return df[(int(df['pos_hg19']) >= start) & (int(df['pos_hg19']) <= end)]
    else:
        # raise error
        print(f"Assembly name {assembly_name} is not supported.")
        return pd.DataFrame()


def get_chromosome_string_for_revel_files_search(chromosome):
    """This function returns the chromosome string for the revel files, that fits to the file names in revel files."""
    if chromosome == 'X' or chromosome == 'Y':
        chromosome_str = f'chrom__{chromosome}'
    elif int(chromosome) < 10:
        chromosome_str = f"chrom_0{chromosome}"
    else:
        chromosome_str = f"chrom_{chromosome}"
    return chromosome_str


def get_relevant_revel_file(path, start, end, chromo_str, assmebly_name):
    # create a list of all paths to the revel data files, that are in the right chromosome.
    chromosome_str = get_chromosome_string_for_revel_files_search(chromosome)
    all_files = os.listdir("C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\benchmarking\\REVEL"
                           "\\revel_prediction_files\\")
    revel_files_chromosome = [file for file in all_files if chromosome_str in file]
    # find the revel file that is in the right start and end locations, there are several files for each chromosome
    revel_file_chromosome_and_location = None

    if assmebly_name == 'GRCh38':
        for file in revel_files_chromosome:
            # Check if the start and end locations are in the column 'grch38_pos' in the file.
            revel_data = pd.read_csv(f"C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\benchmarking\\REVEL"
                                     f"\\revel_prediction_files\\{file}", dtype= {'grch38_pos': str})
            # convert all values in the column to string, to compare them to the start and end locations
            revel_data['grch38_pos'] = revel_data['grch38_pos'].astype(str)
            if start in revel_data['grch38_pos'].values and end in revel_data['grch38_pos'].values:
                revel_file_chromosome_and_location = file
                break
    if assmebly_name == 'GRCh37':
        for file in revel_files_chromosome:
            file_start = int(file.split('_')[-2])
            file_end = int(file.split('_')[-1].split('.')[0])
            if start >= file_start and end <= file_end:
                revel_file_chromosome_and_location = file
                break
    return revel_file_chromosome_and_location


if __name__ == "__main__":

    # Get genes from the dvd file
    dvd_file_path = f"C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\Data\\transmembrane_proteins_data.csv"
    genes = pd.read_csv(dvd_file_path)['gene'].unique().tolist()
    #
    # genes_locations = create_genes_location_dict(genes)
    # 
    # # # save the genes locations to a file
    # # genes_locations_df = pd.DataFrame.from_dict(genes_locations, orient='index')
    # # genes_locations_df.to_csv("C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\benchmarking\\REVEL"
    # #                           "\\genes_locations.csv")

    # Load the genes locations from the file
    genes_locations_df = pd.read_csv("C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\benchmarking\\REVEL"
                                    "\\genes_locations.csv")
    # convert the dataframe to a dictionary
    genes_locations = genes_locations_df.set_index('Unnamed: 0').T.to_dict('list')

    # create revel data files, by cutting the dataframes by the start and end locations
    for gene in genes:
        # Get locations
        chromosome, start, end, strand_num, assmebly_name = genes_locations[gene]
        revel_file = get_relevant_revel_file("C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\benchmarking\\REVEL"
                                                "\\revel_prediction_files\\", start, end, chromosome, assmebly_name)
        if revel_file is None:
            print(f"No REVEL data for {gene}.")
            continue
        else:
            # Cut the dataframe by the start and end locations
            revel_relevant_data = pd.read_csv(f"C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\benchmarking\\REVEL"
                                              f"\\revel_prediction_files\\{revel_file}")
            revel_gene_df = cut_dataframe_by_start_end_locations(revel_relevant_data, start, end)

        if revel_gene_df.empty:
            print(f"No REVEL data for {gene}.")
            continue
        print(revel_gene_df.head())

        # if not revel_file_chromosome_and_location:
        #     print(f"Failed to find the right revel file for {gene}.")
        #     continue
        # else:
        #     revel_data = pd.read_csv(revel_file_chromosome_and_location, header=None)
        #     revel_gene_data = cut_dataframe_by_start_end_locations(revel_data, start, end)
        #     reverse = True if strand_num == -1 else False
        #     # save the new dataframe to csv
        #     if reverse:
        #         revel_gene_data.to_csv(f"{gene}_reverse.csv", index=False)
        #     else:
        #         revel_gene_data.to_csv(f"{gene}.csv", index=False)
        #
        # reverse = True if strand_num == -1 else False
        #
        # if reverse:
        #     revel_gene_data = pd.read_csv(f"{gene}_reverse.csv", header=None)
        # else:
        #     revel_gene_data = pd.read_csv(f"{gene}.csv", header=None)
        #
        # # Add header for the columns
        # revel_data.columns = ['chr', 'pos_hg19', 'pos_grch38', 'ref_na', 'alt_na', 'aa_wt', 'aa_mut', 'REVEL', "transcript_id"]
        #
        # if reverse:
        #     revel_data = revel_data[::-1].reset_index(drop=True)
        #     revel_data = add_aa_position_to_df_reverse(revel_data)
        # else:
        #     revel_data = add_aa_position_to_df(revel_data)
        # revel_data.to_csv(f'{gene}_revel_with_pos.csv', index=False)
