import subprocess as sp
import pandas as pd
import blosum as sm  # blosum matrix, from blosum.py in this directory
import csv
import os
import general_tools as tools
import residue_properties as rp
import uniprot_info as uni

""""
structure of extract_features.py:
# extract_all_features_for_all_variants_in_df
# extract_all_features_for_variant
# write_feature_to_df
# extract_feature_one_value / extract_feature_wt_mut / extract_feature_one_value_for_protein
# get_stability / get_oda / get_opra / get_plddt / get_hydrophobicity / get_aa_volume / get_substitution_matrix_value
# / functions from other modules 

Calling this script from main, to update the features.csv file:
    features_df = pd.read_csv(f"{PATH_TO_OUTPUT_FOLDER}features.csv", header=0)
    ext_feat.main(features_df, PATH_TO_OUTPUT_FOLDER)
"""

FEATURES_LIST = ['stability', 'blosum', 'hydrophobicity', 'volume', 'plddt_residue', 'secondary_structure',
                 'sequence_length', 'protein_contain_transmembrane', 'is_residue_transmembranal', 'aa']

FEATURES_LIST_WITH_DELTA = ['stability', 'hydrophobicity', 'volume']

PATH_TO_ERROR_FILE = "/home/inbar/variants/Benign_for_gene_specific/feature_extraction_errors.csv"


def main(features_df: pd.DataFrame, path_to_output_folder: str):
    """Main function for feature extraction."""
    path_to_csv_file = f"{path_to_output_folder}/Benign_for_gene_specific/benign.csv"
    features_df = extract_all_features_for_all_variants_in_df(features_df, path_to_csv_file)
    features_df.to_csv(path_to_csv_file, index=False)


def extract_all_features_for_all_variants_in_df(features_df: pd.DataFrame, path_to_csv_file: str):
    """Extracts all the features for all the variants in the given dataframe, and adds them to the dataframe.
    Args:
        features_df (pd.DataFrame): The dataframe to add the features to.
    """
    counter = 0

    # If run from the start, get all the variants
    variant_to_add = features_df

    # If run from some point in the middle, get all the variants that don't have features yet
    # variant_to_add = features_df[features_df['stability_WT'].isnull()]

    # # If you want to run a specific variant, get the variant
    # variant = 'N2159L'
    # gene = 'MYO7A'
    # variant_to_add = features_df[(features_df['variant'] == variant) & (features_df['gene'] == gene)]

    for index, row in variant_to_add.iterrows():
        # For re-runs, Check if the variant is already in the dataframe with some features
        gene = row['gene']
        variant = row['variant']
        counter += 1
        # if row['pathogenicity'] == 'benign':
        #     pathogenicity = 'Benign'
        # else:
        #     pathogenicity = 'Pathogenic'
        variant_path = tools.get_variant_path(gene, variant, "Benign_for_gene_specific")
        print(f"{counter}...extracting features for variant {variant} in gene {gene}")

        # Extract all the features for this variant
        # To change the features to extract, add them to the FEATURES_LIST up top.
        # features_df = extract_features_for_variant(variant_path, features_df, gene, variant)

        # Add delta columns using get_delta_feature, for this variant
        # To change/add more features with delta, add them to the FEATURES_LIST_WITH_DELTA list up top.
        features_df = extract_delta_values_for_variant(features_df, gene, variant)

        # Append the variant to the csv
        features_df.to_csv(path_to_csv_file, index=False, header=True)
    return features_df


def extract_features_for_variant(variant_path: str, features_df: pd.DataFrame, gene: str, variant: str):
    """Extracts all the features for the given variant.
    Args:
        variant_path (str): The path to the variant folder.
        features_df (pd.DataFrame): The dataframe to add the features to.
        gene (str): The gene of the variant. Allows to insert the value in the right position in the csv.
        variant (str): The variant. Allows to insert the value in the right position in the csv.
        aa1 (str): The first amino acid in the variant. Allows to insert the value in the right position in the csv.
        aa2 (str): The second amino acid in the variant. Allows to insert the value in the right position in the csv.
    """
    # write all features to dataframe
    try:
        folder_name = tools.get_folder_name_from_path(variant_path)
        aa1 = tools.get_amino_acid_of_wt_from_variant_folder(folder_name)
        aa2 = tools.get_amino_acid_of_variant_from_variant_folder(folder_name)
        for feature in FEATURES_LIST:
            features_df = write_feature_to_df(feature, variant_path, features_df, gene, variant, aa1, aa2)
            print(f"{feature} done!")
    except Exception as e:
        print(f"Error in extract_features_for_variant: {e}")
        # write error to csv with the variant and the error
        write_error_to_csv(gene, variant, e, PATH_TO_ERROR_FILE)
    return features_df


def write_error_to_csv(gene: str, variant: str, error, path_to_errors_csv: str):
    """Writes the error to an errors csv file."""
    # Load the errors csv file, or create it if it doesn't exist
    try:
        errors_df = pd.read_csv(path_to_errors_csv)
    except FileNotFoundError:
        errors_df = pd.DataFrame(columns=['gene', 'variant', 'error'])
    error = str(error)
    errors_df = errors_df.append({'gene': gene, 'variant': variant, 'error': error}, ignore_index=True)
    with open(path_to_errors_csv, 'a') as f:
        errors_df.to_csv(f, header=False)


def get_delta_feature(feature: str, row: pd.Series):
    """Returns the delta value for the given feature, in the given row of the dataframe."""
    return row[f"{feature}_MUT"] - row[f"{feature}_WT"]


def add_delta_columns(df: pd.DataFrame) -> pd.DataFrame:
    """This function adds a column for each two columns with MUT and WT part.
    This column contains the subduction of the value in column MUT from the value in column WT.
    For example,
    df contains "oda_WT" and "oda_MUT".
    This function calculates, for each row, the value oda_MUT - oda_WT,
    and returns this as a column called "oda_delta"."""
    for feature in FEATURES_LIST_WITH_DELTA:
        df[f"{feature}_delta"] = df[f"{feature}_MUT"] - df[f"{feature}_WT"]
    return df


def extract_delta_values_for_variant(features_df: pd.DataFrame, gene: str, variant: str):
    """Extracts the delta values for all the variants in the given dataframe, and adds them to the dataframe."""
    # For every row in the dataframe, add the delta values, using the add_delta_columns function
    for feature in FEATURES_LIST_WITH_DELTA:
        features_df.loc[(features_df['gene'] == gene) & (features_df['variant'] == variant), f"{feature}_delta"] = \
            get_delta_feature(feature, features_df.loc[(features_df['gene'] == gene) & (features_df['variant'] == variant)])
    return features_df


def write_feature_to_df(feature: str, variant_path: str, features_df: pd.DataFrame, gene: str, variant: str,
                        aa1: str = None, aa2: str = None):
    """Adds the given feature to the given dataframe.
    Args:
        feature (str): The feature to extract, e.g. 'blosum', 'hydrophobicity', 'volume'.
        variant_path (str): The path to the variant folder.
        features_df (pd.DataFrame): The dataframe to add the feature to.
        gene (str): The gene of the variant. Allows to insert the value in the right position in the csv.
        variant (str): The variant. Allows to insert the value in the right position in the csv.
        aa1 (str): The first amino acid.
        aa2 (str): The second amino acid.
    """
    print(f"Feature is: {feature} for variant {variant} in gene {gene}")
    if feature in ["blosum", "plddt_residue", "secondary_structure", "sequence_length", "is_residue_transmembranal",
                   "consurf_score"]:
        feature_value = extract_feature_one_value(feature, variant_path, aa1=aa1, aa2=aa2, gene_name=gene)
        # Add the feature to the dataframe, in colume 'feature' and row where gene == gene and variant == variant
        features_df.loc[(features_df['gene'] == gene) & (features_df['variant'] == variant), feature] = feature_value
    elif feature in ["opra", "oda", "sasa", "stability", "hydrophobicity", "volume"]:
        wt_value, mut_value = extract_feature_wt_mut(feature, variant_path, aa1=aa1, aa2=aa2)
        # Add the feature to the dataframe, in colume 'feature' and row where gene == gene and variant == variant
        features_df.loc[(features_df['gene'] == gene) & (features_df['variant'] == variant), feature + "_WT"] = wt_value
        features_df.loc[(features_df['gene'] == gene) & (features_df['variant'] == variant), feature + "_MUT"] = \
            mut_value
    elif feature in ["RSA"]:
        wt_value, mut_value = extract_feature_wt_mut(feature, variant_path, aa1=aa1, aa2=aa2, df=features_df)
        # Add the feature to the dataframe, in colume 'feature' and row where gene == gene and variant == variant
        features_df.loc[(features_df['gene'] == gene) & (features_df['variant'] == variant), feature + "_WT"] = wt_value
        features_df.loc[(features_df['gene'] == gene) & (features_df['variant'] == variant), feature + "_MUT"] = \
            mut_value
    elif feature in ["protein_contain_transmembrane"]:
        feature_value = extract_feature_one_value_for_protein(feature, gene)
        features_df.loc[(features_df['gene'] == gene) & (features_df['variant'] == variant), feature] = feature_value
    elif feature in ["aa"]:
        wt_value, mut_value = aa1, aa2
        features_df.loc[(features_df['gene'] == gene) & (features_df['variant'] == variant), "aa_WT"] = wt_value
        features_df.loc[(features_df['gene'] == gene) & (features_df['variant'] == variant), "aa_MUT"] = \
            mut_value
    return features_df


def extract_feature_one_value_for_protein(feature: str, gene_name: str):
    """Extracts the given feature from the given variant. This is one value for entire protein.
        When there's no need for calculations on the protein itself.
        Supports features: 'protein_contain_transmembrane'.
        Args:
            feature (str): The feature to extract, e.g. 'blosum', 'hydrophobicity', 'volume'.
            gene_name (str): The gene to extract for.
        Returns:
            float: The value of the feature.
    """
    feature_mapping = {
        "protein_contain_transmembrane": uni.is_protein_membranal
    }

    if feature in feature_mapping:
        feature_function = feature_mapping[feature]
        if feature in ['sequence_length', 'protein_contain_transmembrane']:
            return feature_function(gene_name)
    else:
        raise ValueError(f"Feature {feature} is not supported, or not covered by this function."
                         f"This function only supports the following features: {feature_mapping.keys()}")


def extract_feature_one_value(feature: str, variant_path: str, aa1: str = None, aa2: str = None, gene_name: str = None) -> float:
    """Extracts the given feature from the given variant.
    Supports features: 'blosum', 'plddt_residue',
    Args:
        feature (str): The feature to extract, e.g. 'blosum', 'hydrophobicity', 'volume'.
        variant_path (str): The path to the variant folder.
        aa1 (str): The first amino acid.
        aa2 (str): The second amino acid.
        gene_name (str): The gene name.
    Returns:
        float: The value of the feature.
    """
    feature_mapping = {
        'blosum': get_substitution_matrix_value,
        'plddt_residue': get_plddt,
        'secondary_structure': uni.get_secondary_structure,
        'sequence_length': uni.get_sequence_length,
        'is_residue_transmembranal': uni.get_type_of_residue_membrane_or_globular,
        'consurf_score': get_consurf_conservation_score
    }
    if feature in feature_mapping:
        feature_function = feature_mapping[feature]
        if feature in ['plddt_residue']:
            uniprot_id = tools.get_uniprot_id_from_path(variant_path)
            af_file_path = f"{variant_path}/AF_{uniprot_id}.pdb"
            residue_num = tools.get_residue_number_from_variant_folder(tools.get_folder_name_from_path(variant_path))
            return feature_function(residue_num, af_file_path)
        elif feature in ['consurf_score']:
            consurf_path = f"home/inbar/consurf/{gene_name}.txt"
            residue_num = tools.get_residue_number_from_variant_folder(tools.get_folder_name_from_path(variant_path))
            return feature_function(consurf_path, gene_name, residue_num)
        elif feature in ['blosum']:
            return feature_function(aa1, aa2)
        elif feature in ['secondary_structure', 'is_residue_transmembranal']:
            position = tools.get_residue_number_from_variant_folder(tools.get_folder_name_from_path(variant_path))
            return feature_function(gene_name, position)
        elif feature in ['sequence_length']:
            return feature_function(gene_name)
    else:
        raise ValueError(f"Feature {feature} is not supported, or not covered by this function."
                         f"This function only supports the following features: {feature_mapping.keys()}")


def extract_feature_wt_mut(feature: str, variant_path: str, aa1: str = None, aa2: str = None, df: pd.DataFrame = None) -> (float, float):
    """Extracts the given feature from the given variant.
    Supports features: 'opra', 'oda', 'sasa', 'stability', 'hydrophobicity', 'volume', 'RSA'.
    Args:
        feature (str): The feature to extract, e.g. 'blosum', 'hydrophobicity', 'volume'.
        variant_path (str): The path to the variant folder.
        aa1 (str): The first amino acid, for features that need two amino acids.
        aa2 (str): The second amino acid, for features that need two amino acids.
        df (pd.DataFrame): The dataframe with the features, for the RSA feature (to get the asa values).
    Returns:
        tuple of float: The value of the feature.
    """
    feature_mapping = {
        'oda': get_oda,
        'opra': get_opra,
        'sasa': get_sasa,
        'stability': get_stability,
        'hydrophobicity': get_hydrophobicity,
        'volume': get_aa_volume,
        'RSA': get_RSA_of_aa
    }

    if feature in feature_mapping:
        feature_function = feature_mapping[feature]
        if feature in ['hydrophobicity', 'volume']:
            return feature_function(aa1, aa2)
        elif feature in ['oda', 'opra', 'sasa', 'stability']:
            return feature_function(variant_path)
        elif feature in ['RSA']:
            folder_name = tools.get_folder_name_from_path(variant_path)
            gene = tools.get_gene_from_variant_folder(folder_name)
            variant = tools.get_variant_from_variant_folder(folder_name)
            asa_wt = df.loc[(df['gene'] == gene) & (df['variant'] == variant)][['sasa_WT']].values[0][0]
            asa_mut = df.loc[(df['gene'] == gene) & (df['variant'] == variant)][['sasa_MUT']].values[0][0]
            return feature_function(aa1, asa_wt), feature_function(aa2, asa_mut)

    else:
        raise ValueError(f"Feature {feature} is not supported, or not covered by this function."
                         f"This function only supports the following features: {feature_mapping.keys()}")


def get_oda(variant_path: str) -> (float, float):
    """Returns the ODA score for the given variant.
    ODA is optimal docking area, a score that represents potential binding sites, using protein surface desolvation
    energy.
    Args:
        variant_path (str): The path to the variant folder.
    Returns:
        tuple of floats: The ODA score for wt and mut.
    """
    # get all the file names that end with .oda.pdb
    file_names = [file for file in os.listdir(variant_path) if file.endswith(".oda.pdb")]
    residue_number = tools.get_variant_location_from_variant_folder(tools.get_folder_name_from_path(variant_path))
    # convert the files to dataframes
    if len(file_names) != 2:
        print(f"WARNING: There are {len(file_names)} files that end with .oda.pdb in {variant_path}.")
        # write error to lop file:
        with open(os.path.join("/home/inbar/log_files/", 'no_two_oda_files.txt'), 'a') as log_file:
            log_file.write(f"{variant_path}.\n")
        return None, None
    wt_oda_df = pd.DataFrame()
    mut_oda_df = pd.DataFrame()
    for file in file_names:
        if file.startswith('AF'):
            wt_oda_df = tools.convert_pdb_to_dataframe(os.path.join(variant_path, file))
        else:
            mut_oda_df = tools.convert_pdb_to_dataframe(os.path.join(variant_path, file))
    # extract the oda values (in the b factor column) for the given residue number
    wt_oda = wt_oda_df.query(f'residue_number == {residue_number}')['b_factor'].iloc[0]
    mut_oda = mut_oda_df.query(f'residue_number == {residue_number}')['b_factor'].iloc[0]
    return wt_oda, mut_oda


def get_opra(variant_path: str) -> (float, float):
    """Extracts the OPRA score from the given variant.
    OPRA is optimal protein-RNA area calculation. It's useful for identifying RNA binding sites on proteins.
    Values equal to or lower than -1 are considered to be RNA binding sites.
    Args:
        variant_path (str): The path to the variant folder.
    Returns:
        tuple of float: The OPRA scores for wt and mut.
    """
    # get all the file names that end with .opra.pdb
    file_names = [file for file in os.listdir(variant_path) if file.endswith(".opra.pdb")]
    residue_number = tools.get_variant_location_from_variant_folder(tools.get_folder_name_from_path(variant_path))
    # convert the files to dataframes
    if len(file_names) != 2:
        print(f"WARNING: There are {len(file_names)} files that end with .opra.pdb in {variant_path}.")
        # write error to lop file:
        with open(os.path.join("/home/inbar/log_files/", 'no_two_opra_files.txt'), 'a') as log_file:
            log_file.write(f"{variant_path}.\n")
        return None, None
    wt_opra_df = pd.DataFrame()
    mut_opra_df = pd.DataFrame()
    for file in file_names:
        if file.startswith('AF'):
            wt_opra_df = tools.convert_pdb_to_dataframe(os.path.join(variant_path, file))
        else:
            mut_opra_df = tools.convert_pdb_to_dataframe(os.path.join(variant_path, file))
    # extract the opra values (in the b factor column) for the given residue number
    wt_opra = wt_opra_df.query(f'residue_number == {residue_number}')['b_factor'].iloc[0]
    mut_opra = mut_opra_df.query(f'residue_number == {residue_number}')['b_factor'].iloc[0]
    return wt_opra, mut_opra


def get_sasa(variant_path: str) -> (float, float):
    """Gets the path for a variant folder, and returns the SASA of the WT and MUT proteins.
    Uses both files created by run_sasa.
    The sasa value is the sum of the sasa values of all the atoms in the residue.
    Returns a tuple of (WT SASA, MUT SASA)"""
    file_names = [file for file in os.listdir(variant_path) if file.endswith(".asa.pdb")]
    residue_number = tools.get_variant_location_from_variant_folder(tools.get_folder_name_from_path(variant_path))
    # convert the files to dataframes
    if len(file_names) != 2:
        print(f"WARNING: There are {len(file_names)} files that end with .asa.pdb in {variant_path}.")
        # write error to lop file:
        with open(os.path.join("/home/inbar/log_files/", 'no_two_sasa_files.txt'), 'a') as log_file:
            log_file.write(f"{variant_path}.\n")
        return None, None
    wt_sasa_df = pd.DataFrame()
    mut_sasa_df = pd.DataFrame()
    for file in file_names:
        if file.startswith('AF'):
            wt_sasa_df = tools.convert_pdb_to_dataframe(os.path.join(variant_path, file))
        else:
            mut_sasa_df = tools.convert_pdb_to_dataframe(os.path.join(variant_path, file))
    # extract the sasa values (in the b factor column) for the given residue number
    wt_sasa_for_residue = wt_sasa_df.query(f'residue_number == {residue_number}')['b_factor']
    mut_oda_for_residue = mut_sasa_df.query(f'residue_number == {residue_number}')['b_factor']
    # sum the sasa values
    wt_sasa = sum(wt_sasa_for_residue)
    mut_sasa = sum(mut_oda_for_residue)
    return wt_sasa, mut_sasa


def get_stability(folder_path: str) -> (float, float, float):
    """Gets the path for a variant folder, and returns the stability of the WT and MUT proteins.
    This is extracted from the text file that starts with the word energies in the folder.
    For example, file name could be: energies_341_AF_Q9HCJ1.txt

    Returns a tuple of (WT stability, MUT stability)"""
    # get the file name, which starts with the word energies - there is only one such file in the folder
    file_name = [file for file in os.listdir(folder_path) if file.startswith("energies")][0]
    # get the stability
    with open(os.path.join(folder_path, file_name), 'r') as f:
        lines = f.readlines()
    # get the stability
    wt_stability = float(lines[0].split()[1])
    mut_stability = float(lines[1].split()[1])
    return wt_stability, mut_stability


def get_hydrophobicity(aa: str, aa2: str) -> (float, float):
    """Returns the hydrophobicity of the residue"""
    return rp.get_aa_hydrophobicity_by_kd_scale(aa), rp.get_aa_hydrophobicity_by_kd_scale(aa2)


def get_aa_volume(aa: str, aa2: str) -> (float, float):
    """Returns the volume of the residue"""
    return rp.get_aa_volume(aa), rp.get_aa_volume(aa2)


def get_hydrogen_bonds():
    """Returns the number of hydrogen bonds of residue"""
    pass


def convert_pdb_to_csv(file_name: str):
    """Converts the given PDB file to a CSV file, and returns the path to the CSV file."""
    # Check if the file is a PDB file
    if not file_name.endswith(".pdb"):
        raise ValueError("The given file is not a PDB file.")

    # read the pdb file
    with open(file_name, 'r') as pdb_file:
        lines = pdb_file.readlines()
    # create a new csv file
    with open(f"{file_name}.csv", 'w') as csv_file:
        writer = csv.writer(csv_file)
        for line in lines:
            if line.startswith("ATOM"):
                writer.writerow(line.split())

    return f"{file_name}.csv"


def get_plddt(residue_num: str, af_file_path: str):
    """returns the pLDDT score of the residue
    Uses the alphafold pdb file to extract the pLDDT score"""
    residue_num = int(residue_num)
    # convert file to csv
    af_df = tools.convert_pdb_to_dataframe(af_file_path)
    # get the pLDDT score, which is the first value in the b_factor column where the residue number = residue_num
    # in the dataframe
    matching_rows = af_df.loc[af_df['residue_number'] == residue_num, 'b_factor']
    if not matching_rows.empty:
        value = matching_rows.iloc[0]
        print(value)
        return value
    else:
        print(f"No pLDDT score found for residue number {residue_num}.")
        return None


def get_average_plddt_score(pdb_file_pwd: str):
    pass


def get_substitution_matrix_value(wt_aa: str, mut_aa: str) -> int:
    """Returns the substitution value of the given mutation.
    Uses the BLOSUM62 matrix.
    Args:
        wt_aa: the wild type amino acid
        mut_aa: the mutated amino acid
    Returns:
        int: the substitution value of the given mutation, according to the BLOSUM62 matrix.
    """
    return sm.get_blosum62_value(wt_aa, mut_aa)


# def get_consurf_conservation_score_for_all_residues(path_to_gene_folder: str) -> dict:
#     """Returns the conservation score of all the residues in the protein.
#     Uses the consurf file, which is a text file downloaded from the consurf website.
#     The conservation score is the 5th column in the file, and the residue number is in the first column.
#     The function finds the value in the 5th column and in the row that corresponds to the given residue number.
#     It skips the beginning of the file, which ends with a line that starts with "POS".
#     Args:
#         path_to_gene_folder: the path to the gene folder, to get the consurf file
#     Returns:
#         dict: a dictionary of {residue number: conservation score}
#     """
#     # get the consurf file
#     consurf_file = [file for file in os.listdir(path_to_gene_folder) if file.startswith("consurf")][0]
#     # read the file
#     with open(os.path.join(path_to_gene_folder, consurf_file), 'r') as f:
#         lines = f.readlines()
#     # get the conservation score for each residue
#     conservation_scores = {}
#     for line in lines:
#         if line.startswith(" POS") or line.startswith("  "):
#             print("line ", line)
#             continue
#     for line in lines:
#         line_parts = line.split()
#         print("line parts ", line_parts)
#         if len(line_parts) > 5:
#             conservation_scores[line_parts[0]] = line_parts[4]
#     print(conservation_scores)
#     return conservation_scores


def get_conserf_coservation_scores_for_gene(path_to_consurf_folder: str, gene_name: str) -> dict:
    """Return the conservation scores of all the residues in the gene.
    Searches
    Uses the consurf file, which is a text file downloaded from the consurf website.
    In the path, the consurf folder should contain a file with the name "gene_name.txt".
    The conservation score is the 5th column in the file, and the residue number is in the first column.
    The function finds the value in the 5th column and in the row that corresponds to the given residue number.
    It skips the beginning of the file, which ends with a line that starts with "POS".
    Args:
        path_to_consurf_folder: the path to the consurf folder, which should contain a txt file with the name
            "gene_name.txt"
        gene_name: the name of the gene
    Returns:
        dict: a dictionary of {residue number: conservation score}
    """
    # check if the path is valid
    if not os.path.isdir(path_to_consurf_folder):
        raise ValueError("The given path is not a valid path to a folder.")
    # Check if the folder contains a file with the name "gene_name.txt"
    if not any([file.startswith(gene_name) for file in os.listdir(path_to_consurf_folder)]):
        raise ValueError(f"The given folder does not contain a file with the name {gene_name}.")

    # get the consurf file with the name "gene_name.txt"
    consurf_file = [file for file in os.listdir(path_to_consurf_folder) if file.startswith(gene_name)][0]
    # read the file
    with open(os.path.join(path_to_consurf_folder, consurf_file), 'r') as f:
        lines = f.readlines()
    # get the conservation score for each residue
    conservation_scores = {}
    for line in lines:
        if line.startswith(" POS") or line.startswith("  "):
            print("line ", line)
            continue
    for line in lines:
        line_parts = line.split()
        print("line parts ", line_parts)
        if len(line_parts) > 5:
            conservation_scores[line_parts[0]] = line_parts[4]
    print(conservation_scores)
    return conservation_scores


def get_consurf_conservation_score_for_residue(conservation_scores_dict: dict, residue_num: str) -> int:
    """Returns the conservation score of the given residue.
    Uses the consurf file, which is a text file downloaded from the consurf website.
    The conservation score is the 5th column in the file, and the residue number is in the first column.
    The function finds the value in the 5th column and in the row that corresponds to the given residue number.
    It skips the beginning of the file, which ends with a line that starts with "POS".
    Args:
        conservation_scores_dict: a dictionary of {residue number: conservation score}
        residue_num: the residue number
    Returns:
        int: the conservation score of the given residue
    """
    # Check if the residue number is in the dictionary
    if residue_num not in conservation_scores_dict.keys():
        print(f"No conservation score found for residue number {residue_num}.")
        raise ValueError(f"No conservation score found for residue number {residue_num}.")
    return conservation_scores_dict[residue_num]


def get_consurf_conservation_score(consurf_folder: str, gene_name: str, residue_num: str) -> int:
    """Returns the conservation score of the given residue.
    Uses the consurf file, which is a text file downloaded from the consurf website.
    The conservation score is the 5th column in the file, and the residue number is in the first column.
    The function finds the value in the 5th column and in the row that corresponds to the given residue number.
    It skips the beginning of the file, which ends with a line that starts with "POS".
    Args:
        consurf_folder: the path to the consurf folder, which should contain a txt file with the name "gene_name.txt"
        gene_name: the name of the gene
        residue_num: the residue number
    Returns:
        int: the conservation score of the given residue
    """
    # get the conservation scores dictionary
    conservation_scores_dict = get_conserf_coservation_scores_for_gene(consurf_folder, gene_name)
    # get the conservation score for the given residue
    return get_consurf_conservation_score_for_residue(conservation_scores_dict, residue_num)


# def get_consurf_conservation_score(path_to_gene_folder: str, path_to_variant_folder: str) -> int:
#     """Returns the conservation score of the given residue.
#     Uses the consurf file, which is a text file downloaded from the consurf website.
#     The conservation score is the 5th column in the file, and the residue number is in the first column.
#     The function finds the value in the 5th column and in the row that corresponds to the given residue number.
#     It skips the beginning of the file, which ends with a line that starts with "POS".
#     Args:
#         path_to_gene_folder: the path to the gene folder, to get the consurf file
#         path_to_variant_folder: the path to the variant folder, to get the residue number
#     Returns:
#         float: the conservation score of the given residue
#
#     Use:
#         consurf = ext_feat.get_consurf_conservation_score("/home/inbar/check/ACTB/", "/home/inbar/check/ACTB/ACTB_P60709_Q189R")
#     """
#     # get the folder name from the path
#     folder_name = tools.get_folder_name_from_path(path_to_variant_folder)
#     # get the residue number from the folder name
#     residue_num = tools.get_residue_number_from_variant_folder(folder_name)
#     residue_num = str(residue_num)  # convert to string
#     score = 0  # the score to be returned
#     # read the consurf file
#     consurf_file_path = os.path.join(path_to_gene_folder, "consurf_grades.txt")
#     with open(consurf_file_path, 'r') as f:
#         for line in f:
#             if line.startswith(" POS"):
#                 break
#         for line in f:
#             # find the line that starts with the residue number
#             if line.strip().startswith(residue_num):
#                 # get the conservation score
#                 score = line.split()[4]
#                 break
#
#     # if the residue number was not found, write an error message to the log file
#     if score == 0:
#         with open("log.txt", 'a') as log_file:
#             log_file.write(f"Residue number {residue_num} in path {path_to_gene_folder} was not found in the consurf"
#                            f" file.\n")
#     if "*" in str(score):  # if the score has a star, it means that the conservation score is not reliable
#         return 10  # 10 stands for not reliable
#     else:
#         return int(score)


# Feature engineering #

def get_RSA_of_aa(aa: str, asa: float) -> float:
    """Returns the relative surface area (RSA) of the given amino acid.
    The RSA is calculated by dividing the absolute surface area (ASA) by the maximum ASA of the given amino acid.
    According to the paper "Maximum Allowed Solvent Accessibilites of Residues in Proteins" by Tien et al., 2013.
    RSA = ASA / max ASA
    Args:
        aa: the amino acid
        asa: the absolute surface area of the amino acid
    Returns:
        float: the RSA of the given amino acid
    """
    # Dictionary for the relative surface area of each amino acid
    max_asa = {'A': 121, 'R': 265, 'N': 187, 'D': 187, 'C': 148, 'Q': 214, 'E': 214, 'G': 97, 'H': 216, 'I': 195,
               'L': 191, 'K': 230, 'M': 203, 'F': 228, 'P': 154, 'S': 143, 'T': 163, 'W': 264, 'Y': 255, 'V': 165}
    rsa = asa / max_asa[aa]
    return rsa

