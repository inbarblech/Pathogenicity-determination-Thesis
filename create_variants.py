import uniprot_info as uni
import pandas as pd

TAGS = {1: 'pathogenic', 0: 'benign'}

def create_variants_df_from_arff(path_to_arff: str) -> pd.DataFrame:
    """Create a dataframe from the given arff file."""


def create_benign_variants_df_from_arff(path_to_arff: str, path_to_benign_csv: str) -> pd.DataFrame:
    """Create a dataframe from the given arff file.
    The dataframe will contain only the benign variants, from the PatMut file."""
    list_of_genes = ['GJB2']

    uniprot_of_genes = {gene: uni.get_uniprot_id(gene) for gene in list_of_genes}
    list_of_uniprot_ids = list(uniprot_of_genes.values())

    dataframes = []
    for uniprot_id in list_of_uniprot_ids:
        path_to_arff = f"C:\\Users\\InbarBlech\\Downloads\\{uniprot_id}.arff"
        df = pd.read_csv(path_to_arff, skiprows=13, header=None, sep=',', engine='python')
        df.columns = ['variant', 'uniprot_id', 'vdw_volume', 'hydrophobicity', 'substitution_matrix', 'pssm_native',
                      'entropy', 'imp_res', 'tag']
        # Add gene column, by using the uniprot_of_genes dictionary
        benign_df = df[df['tag'] == '0']
        # change dataframe to fit the format of the other dataframes: gene, variant, pathogenicity, uniprot_id
        benign_df['gene'] = benign_df['uniprot_id'].map({v: k for k, v in uniprot_of_genes.items()})
        # write "benign" in the pathogenicity column
        benign_df.loc[:, 'pathogenicity'] = 'benign'
        benign_df = benign_df[['gene', 'variant', 'pathogenicity']]
        benign_df = benign_df.drop_duplicates()
        benign_df = benign_df.reset_index(drop=True)
        # add uniprot id column
        benign_df['uniprot_id'] = uniprot_id
        # add to list of dataframes
        dataframes.append(benign_df)

    # create one dataframe with all the benign variants
    benign_df = pd.DataFrame(columns=['gene', 'variant', 'pathogenicity', 'uniprot_id'])
    for dataframe in dataframes:
        benign_df = benign_df.append(dataframe, ignore_index=True)

    # write to file
    benign_df.to_csv(f"C:\\Users\\InbarBlech\\OneDrive - mail.tau.ac.il\\Documents\\\Thesis\\Classification project\\"
                     f"Data\\benign_artificial_variants_patmut\\benign.csv", index=False)
    return benign_df


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

    create_benign_variants_df_from_arff("C:\\Users\\InbarBlech\\Downloads\\P29033.arff")

