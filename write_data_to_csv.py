import pandas as pd
import os
import dvd_data as dvd
import uniprot_info as uni
import dvd_data as dvd

PATH_TO_DATA_FOLDER = "/home/inbar/DVDdata/"
PATH_TO_CSV_BENIGN = "/home/inbar/all_BENIGN_variants.csv"
PATH_TO_CSV_PATHOGENIC = "/home/inbar/all_PATHOGENIC_variants.csv"
PATH_TO_VARIANTS_FOLDER = "/home/inbar/variants/"
PATHOGENICITY_TYPES = ["Benign", "Pathogenic"]


def create_new_df():
    df = pd.DataFrame(columns=['gene', 'pathogenicity', 'uniprot_id', 'variant', 'pathogenicity', 'delta_stability',
                               'delta_energy'])
    return df


def add_variant_data_to_df(df, variant_data: dict):
    """Adds the variant data to the csv file."""
    df = df.append(variant_data, ignore_index=True)
    return df


def add_feature_to_csv(variant, gene, feature_value, feature, csv_path) -> None:
    """Add feature to csv file.
    Save the feature value in the csv file at the variant and gene location.
    Does not create a new csv file, but rather updates the existing one.
    Args:
        variant (str): The variant name.
        gene (str): The gene name.
        feature_value (float): The value of the feature.
        feature (str): The name of the feature.
        csv_path (str): The path to the csv file.
    """
    df = pd.read_csv(csv_path)
    # add feature to df at variant and gene location:
    df.loc[(df['Variant'] == variant) & (df['Gene_name'] == gene), feature] = feature_value
    df.to_csv(csv_path, index=False)


# This function created the csv file with all the variants, not needed after the csv file is created.
def create_variants_csv_from_folder(path_to_folder_pathogenicity, output_filename):
    """Creates a csv file with all the variants in the given folder."""
    df = pd.DataFrame(columns=['gene', 'pathogenicity', 'variant', 'uniprot_id'])
    pathogenicity = os.path.basename(path_to_folder_pathogenicity)
    for root, dirs, files in os.walk(path_to_folder_pathogenicity):
        counter = 0
        for file_name in files:
            if file_name.endswith(".csv"):
                file_path = os.path.join(root, file_name)
                gene_name = dvd.extract_gene_name(file_name)
                uniprot_id = uni.get_uniprot_id(gene_name)
                variants = dvd.make_list_of_variants_per_gene_from_dvd_data(file_path)
                for variant in variants:
                    df = df.append({'gene': gene_name, 'pathogenicity': pathogenicity, 'variant': variant,
                                    'uniprot_id': uniprot_id}, ignore_index=True)
                    print(f"{counter} - Variant {variant} from gene {gene_name} added to DataFrame.")
    output_path = os.path.join(os.getcwd(), output_filename)  # Construct the output path using os.getcwd()
    df.drop_duplicates(inplace=True)  # Remove duplicate rows from the DataFrame (if they exist)
    df.to_csv(output_path, index=False)  # Save the DataFrame as CSV without including the index
    print(f"CSV file saved to: {output_path}")


"""Some preprocessing of the raw data:"""


def remove_duplicate_rows_from_csv(csv_path):
    """Removes identical rows from a csv file, if they exist."""
    df = pd.read_csv(csv_path)
    df.drop_duplicates(inplace=True)
    df.to_csv(csv_path, index=False)


def find_duplicate_rows_in_csv(csv_path):
    """Finds duplicate rows in a csv file."""
    df = pd.read_csv(PATH_TO_CSV_PATHOGENIC)
    duplicate_rows = df[df.duplicated()]
    print(duplicate_rows)
    print(len(duplicate_rows))