import os

from data_retrievel_and_feature_extraction import uniprot_info as uni
import pandas as pd

TAGS = {1: 'pathogenic', 0: 'benign'}

def create_variants_df_from_arff(path_to_arff: str) -> pd.DataFrame:
    """Create a dataframe from the given arff file."""


def create_benign_variants_df_from_arff(uniprot_ids, path_to_arff: str, path_to_benign_csv: str) -> pd.DataFrame:
    """Create a dataframe from the given arff file.
    The dataframe will contain only the benign variants, from the PatMut file."""

    with open("C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\uniprot_and_genes.txt", "r") as file:
        # The file is in the format: gene, uniprot_id
        uniprot_of_genes = {line.split(",")[0]: line.split(",")[1].strip() for line in file}

    dataframes = pd.DataFrame(columns=['variant', 'pathogenicity', 'uniprot_id'])
    for uniprot_id in uniprot_ids:
        # Check if path to arff file exists
        if not os.path.exists(f"{path_to_arff}\\{uniprot_id}.arff"):
            print(f"Path to arff file does not exist: {path_to_arff}, gene: {uniprot_id}")
            continue
        else:
            df = pd.read_csv(f"{path_to_arff}\\{uniprot_id}.arff", skiprows=13, header=None, sep=',', engine='python')
            df.columns = ['variant', 'uniprot_id', 'vdw_volume', 'hydrophobicity', 'substitution_matrix', 'pssm_native',
                          'entropy', 'imp_res', 'tag']
            # Create a dataframe with only the benign variants, from the PatMut file
            benign_df = df[df['tag'] == 0]
            # Remove all columns except for the uniprot_id, variant, and tag columns
            benign_df = benign_df[['uniprot_id', 'variant', 'tag']]
            # Change the tag column name to pathogenicity
            benign_df = benign_df.rename(columns={'tag': 'pathogenicity'})
            # Remove duplicates
            benign_df = benign_df.drop_duplicates()
            benign_df = benign_df.reset_index(drop=True)
            # Add to dataframes dataframe
            dataframes = pd.concat([dataframes, benign_df], ignore_index=True)
            print(f"Finished gene: {uniprot_id}")

    # write to file
    dataframes.to_csv(path_to_benign_csv, index=False)
    return dataframes


def get_dataframe_without_overlap(current_variants: pd.DataFrame, new_variants: pd.DataFrame) -> pd.DataFrame:
    """Function recieves two dataframes: current_variants and new_variants.
    It checks if there are variants in new_variants that already exist in current_variants.
    If there are, it prints them to the console, and returns a new dataframe with only the new variants.
    If there aren't, it returns the new_variants dataframe as is.

    It checks for overlap by comparing the variant and gene columns.
    """
    # check if there are variants in new_variants that already exist in current_variants
    # if there are, print them to the console
    overlap = new_variants[new_variants[['variant', 'gene']].isin(current_variants[['variant', 'gene']]).all(axis=1)]
    if not overlap.empty:
        print(f"Overlap between current variants and new variants: {overlap}")
        # return a new dataframe with only the new variants
        new_variants = new_variants[~new_variants[['variant', 'gene']].isin(current_variants[['variant', 'gene']]).all(axis=1)]
        new_variants = new_variants.reset_index(drop=True)
        return new_variants
    else:
        return new_variants


def remove_duplicates_from_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """Function recieves a dataframe, and returns a new dataframe with no duplicates."""
    df = df.drop_duplicates()
    df = df.reset_index(drop=True)
    return df


if __name__ == "__main__":
    # new_var = pd.read_csv("C:\\Users\\InbarBlech\\OneDrive - mail.tau.ac.il\\Documents\\Thesis"
    #                       "\\Classification project\\Data\\benign_artificial_variants_patmut\\GJB2_benign.csv")
    # existing_var = pd.read_csv("C:\\Users\\InbarBlech\\OneDrive - mail.tau.ac.il\\Documents\\Thesis\\Findings\\"
    #                            "features.csv")
    # print(f"First, length is {len(new_var)}")
    #
    # new_var = get_dataframe_without_overlap(existing_var, new_var)
    # print(f"After reducing duplicates, length is {len(new_var)}")
    #
    # # write to file
    # new_var.to_csv("C:\\Users\\InbarBlech\\OneDrive - mail.tau.ac.il\\Documents\\Thesis\\Classification project\\Data\\"
    #                  "\\benign_artificial_variants_patmut\\benign_without_duplicates_with_features.csv", index=False)

    # Load data from all dvd
    dvd_df = pd.read_csv("C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\Data\\data_for_all_dvd.csv")
    # Add uniprot_id column to the dataframe, by using the uniprot_and_genes dictionary
    # read the dictionary from the file
    with open("C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\uniprot_and_genes.txt", "r") as file:
        # The file is in the format: gene, uniprot_id
        uniprot_of_genes = {line.split(",")[0]: line.split(",")[1].strip() for line in file}
    dvd_df['uniprot_id'] = dvd_df['gene'].map(uniprot_of_genes)

    # Save the dataframe to a csv file
    dvd_df.to_csv("C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\Data\\data_for_all_dvd_with_uniprot_id.csv", index=False)

    # Create a list of Uniprot IDs
    list_of_uniprot_ids = dvd_df['uniprot_id'].tolist()
    # Remove duplicates
    list_of_uniprot_ids = list(dict.fromkeys(list_of_uniprot_ids))
    print(f"Number of uniprot ids: {len(list_of_uniprot_ids)}")
    # Create benign variants df
    benign_df = pd.DataFrame(columns=['gene', 'variant', 'pathogenicity', 'uniprot_id'])
    benign_df = create_benign_variants_df_from_arff(list_of_uniprot_ids , f"C:\\Users\\InbarBlech\\Downloads\\"
                                                            f"patmut_all\\patmut_20230905",
                                                            "C:\\Users\\InbarBlech\\PycharmProjects"
                                                            "\\Thesis\\Data\\patmut\\benign_variants_patmut_25_4_24.csv")

