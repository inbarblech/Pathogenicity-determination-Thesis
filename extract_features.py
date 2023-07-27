import subprocess as sp
import pandas as pd
import blosum as sm  # blosum matrix, from blosum.py in this directory
import csv
import os
import general_tools as tools
import residue_properties as rp

""""
structure of extract_features.py:
# extract_all_features_for_all_variants_in_df
# extract_all_features_for_variant
# write_feature_to_df
# extract_feature_one_value / extract_feature_wt_mut
# get_stability / get_oda / get_opra / get_plddt / get_hydrophobicity / get_aa_volume / get_substitution_matrix_value
"""


def extract_all_features_for_all_variants_in_df(features_df: pd.DataFrame):
    """Extracts all the features for all the variants in the given dataframe, and adds them to the dataframe.
    Args:
        features_df (pd.DataFrame): The dataframe to add the features to.
    """
    for index, row in features_df.iterrows():
        gene = row['gene']
        variant = row['variant']
        pathogenicity = row['pathogenicity']
        variant_path = tools.get_variant_path(gene, variant, pathogenicity)
        print(f"extracting features for variant {variant} in gene {gene} in pathogenicity {pathogenicity}")
        extract_all_features_for_variant(variant_path, features_df, gene, variant)


def extract_all_features_for_variant(variant_path: str, features_df: pd.DataFrame, gene: str, variant: str):
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
    folder_name = tools.get_folder_name_from_path(variant_path)
    aa1 = tools.get_amino_acid_of_wt_from_variant_folder(folder_name)
    aa2 = tools.get_amino_acid_of_variant_from_variant_folder(folder_name)
    write_feature_to_df('blosum', variant_path, features_df, gene, variant, aa1=aa1, aa2=aa2)
    print("blosum done")
    write_feature_to_df('hydrophobicity', variant_path, features_df, gene, variant, aa1=aa1, aa2=aa2)
    print("hydrophobicity done")
    write_feature_to_df('volume', variant_path, features_df, gene, variant, aa1=aa1, aa2=aa2)
    print("volume done")
    write_feature_to_df('plddt', variant_path, features_df, gene, variant)
    print("plddt done")
    # write_feature_to_df('opra', variant_path, features_df, gene, variant)
    # write_feature_to_df('oda', variant_path, features_df, gene, variant)
    # write_feature_to_df('sasa', variant_path, features_df, gene, variant)


def write_feature_to_df(feature: str, variant_path: str, features_df: pd.DataFrame, gene: str, variant: str, wt_or_mut: str = None,
                        aa1: str = None, aa2: str = None):
    """Adds the given feature to the given dataframe.
    Args:
        feature (str): The feature to extract, e.g. 'blosum', 'hydrophobicity', 'volume'.
        variant_path (str): The path to the variant folder.
        features_df (pd.DataFrame): The dataframe to add the feature to.
        gene (str): The gene of the variant. Allows to insert the value in the right position in the csv.
        variant (str): The variant. Allows to insert the value in the right position in the csv.
        wt_or_mut (str): 'wt' or 'mut'.
        aa1 (str): The first amino acid.
        aa2 (str): The second amino acid.
    """
    if feature in ["blosum", "plddt"]:
        feature_value = extract_feature_one_value(feature, aa1, aa2)
        # Add the feature to the dataframe, in colume 'feature' and row where gene == gene and variant == variant
        features_df.loc[(features_df['gene'] == gene) & (features_df['variant'] == variant), feature] = feature_value
    elif feature in ["opra", "oda", "sasa", "stability", "hydrophobicity", "volume"]:
        wt_value, mut_value = extract_feature_wt_mut(feature, variant_path)
        # Add the feature to the dataframe, in colume 'feature' and row where gene == gene and variant == variant
        features_df.loc[(features_df['gene'] == gene) & (features_df['variant'] == variant), feature + "_WT"] = wt_value
        features_df.loc[(features_df['gene'] == gene) & (features_df['variant'] == variant), feature + "_MUT"] = \
            mut_value


def extract_feature_one_value(feature: str, aa1: str = None, aa2: str = None) -> float:
    """Extracts the given feature from the given variant.
    Supports features: 'blosum', 'plddt_residue'.
    Args:
        feature (str): The feature to extract, e.g. 'blosum', 'hydrophobicity', 'volume'.
        variant_path (str): The path to the variant folder.
        aa1 (str): The first amino acid.
        aa2 (str): The second amino acid.
    Returns:
        float: The value of the feature.
    """
    feature_mapping = {
        'blosum': get_substitution_matrix_value,
        'plddt_residue': get_plddt,
    }

    if feature in feature_mapping:
        feature_function = feature_mapping[feature]
        if feature in ['plddt_residue']:
            return feature_function(aa1)
        elif feature in ['blosum']:
            return feature_function(aa1, aa2)
    else:
        raise ValueError(f"Feature {feature} is not supported, or not covered by this function."
                         f"This function only supports the following features: {feature_mapping.keys()}")


def extract_feature_wt_mut(feature: str, variant_path: str, aa1: str = None, aa2: str = None) -> (float, float):
    """Extracts the given feature from the given variant.
    Supports features: 'opra', 'oda', 'sasa', 'stability', 'hydrophobicity', 'volume'.
    Args:
        feature (str): The feature to extract, e.g. 'blosum', 'hydrophobicity', 'volume'.
        variant_path (str): The path to the variant folder.
        aa1 (str): The first amino acid, for features that need two amino acids.
        aa2 (str): The second amino acid, for features that need two amino acids.
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
    }

    if feature in feature_mapping:
        feature_function = feature_mapping[feature]
        if feature in ['hydrophobicity', 'volume']:
            return feature_function(variant_path, aa1, aa2)
        elif feature in ['oda', 'opra', 'sasa', 'stability']:
            return feature_function(variant_path)
    else:
        raise ValueError(f"Feature {feature} is not supported, or not covered by this function."
                         f"This function only supports the following features: {feature_mapping.keys()}")


def run_oda(variant_path: str):
    """Creates a file with the ODA score for the given PDB file.
    The file is created in the variant folder.
    Args:
        variant_path (str): The path to the variant folder.
    """
    # change the working directory to the variant folder
    os.chdir(variant_path)

    # get gene_id from folder name
    folder_name = tools.get_folder_name_from_path(variant_path)
    variant_location = tools.get_variant_location_from_variant_folder(folder_name)
    gene_id = folder_name.split('_')[1]

    # Run oda for wt and mut
    amino_acid = tools.get_amino_acid_of_variant_from_variant_folder(folder_name)
    three_letter_amino_acid = tools.convert_1_letter_aa_to_3_letter(amino_acid)
    command = f"pyDock3 {three_letter_amino_acid}{variant_location}_AF_{gene_id}.pdb oda"
    sp.run(command, shell=True)
    command = f"pyDock3 AF_{gene_id}.pdb oda"
    sp.run(command, shell=True)


def get_oda(variant_path: str) -> (float, float):
    """Returns the ODA score for the given variant.
    ODA is optimal docking area, a score that represents potential binding sites, using protein surface desolvation
    energy.
    Args:
        variant_path (str): The path to the variant folder.
    Returns:
        tuple of floats: The ODA score for wt and mut.
    """
    # get all the file names that end with .pdb.oda
    file_names = [file for file in os.listdir(variant_path) if file.endswith(".pdb.oda")]
    residue_number = tools.get_variant_location_from_variant_folder(tools.get_folder_name_from_path(variant_path))
    # convert the files to dataframes
    if len(file_names) != 2:
        print(f"WARNING: There are {len(file_names)} files that end with .pdb.oda in {variant_path}.")
        return None, None
    for file in file_names:
        if file.startswith('AF'):
            wt_oda_df = tools.convert_pdb_to_dataframe(os.path.join(variant_path, file))
        else:
            mut_oda_df = tools.convert_pdb_to_dataframe(os.path.join(variant_path, file))
    # extract the oda values (in the b factor column) for the given residue number
    wt_oda = wt_oda_df.loc[wt_oda_df['residue_number'] == residue_number, 'b_factor'].values[0]
    mut_oda = mut_oda_df.loc[mut_oda_df['residue_number'] == residue_number, 'b_factor'].values[0]
    return wt_oda, mut_oda


def run_opra(variant_path: str, wt_or_mut: str):
    """creates a file with the OPRA score for the given PDB file."""
    # change the working directory to the variant folder
    os.chdir(variant_path)

    # get gene_id from folder name
    folder_name = tools.get_folder_name_from_path(variant_path)
    variant_location = tools.get_variant_location_from_variant_folder(folder_name)
    gene_id = folder_name.split('_')[1]

    if wt_or_mut == 'mut':
        amino_acid = tools.get_amino_acid_of_variant_from_variant_folder(folder_name)
        three_letter_amino_acid = tools.convert_1_letter_aa_to_3_letter(amino_acid)
        command = f"pyDock3 {three_letter_amino_acid}{variant_location}_AF_{gene_id}.pdb opra"
    elif wt_or_mut == 'wt':
        command = f"pyDock3 AF_{gene_id}.pdb opra"
    else:
        raise ValueError("wt_or_mut must be 'wt' or 'mut'")
    sp.run(command, shell=True)


def get_opra(variant_path: str) -> (float, float):
    """Extracts the OPRA score from the given variant.
    OPRA is optimal protein-RNA area calculation. It's useful for identifying RNA binding sites on proteins.
    Values equal to or lower than -1 are considered to be RNA binding sites.
    Args:
        variant_path (str): The path to the variant folder.
    Returns:
        tuple of float: The OPRA scores for wt and mut.
    """
    # get all the file names that end with .pdb.opra
    file_names = [file for file in os.listdir(variant_path) if file.endswith(".pdb.opra")]
    residue_number = tools.get_variant_location_from_variant_folder(tools.get_folder_name_from_path(variant_path))
    # convert the files to dataframes
    if len(file_names) != 2:
        print(f"WARNING: There are {len(file_names)} files that end with .pdb.opra in {variant_path}.")
        return None, None
    for file in file_names:
        if file.startswith('AF'):
            wt_opra_df = tools.convert_pdb_to_dataframe(os.path.join(variant_path, file))
        else:
            mut_opra_df = tools.convert_pdb_to_dataframe(os.path.join(variant_path, file))
    # extract the opra values (in the b factor column) for the given residue number
    wt_opra = wt_opra_df.loc[wt_opra_df['residue_number'] == residue_number, 'b_factor'].values[0]
    mut_opra = mut_opra_df.loc[mut_opra_df['residue_number'] == residue_number, 'b_factor'].values[0]
    return wt_opra, mut_opra


def run_sasa(path_to_variant_folder: str, wt_or_mut: str):
    """Creates a file with the SASA score for the given PDB file."""
    os.chdir(path_to_variant_folder)
    # get gene_id from folder name
    folder_name = tools.get_folder_name_from_path(path_to_variant_folder)
    variant_location = tools.get_variant_location_from_variant_folder(folder_name)
    gene_id = folder_name.split('_')[1]

    if wt_or_mut == 'wt':
        command = f"dr_sasa -m 0 -i AF_{gene_id}.pdb"
    elif wt_or_mut == 'mut':
        amino_acid = tools.get_amino_acid_of_variant_from_variant_folder(folder_name)
        three_letter_amino_acid = tools.convert_1_letter_aa_to_3_letter(amino_acid)
        command = f"dr_sasa -m 0 -i {three_letter_amino_acid}{variant_location}_AF_{gene_id}.pdb"
    else:
        raise ValueError("wt_or_mut must be 'wt' or 'mut'")
    print(command)
    sp.run(command, shell=True)


def get_sasa(variant_path: str) -> (float, float):
    """Gets the path for a variant folder, and returns the SASA of the WT and MUT proteins.
    Uses both files created by run_sasa.
    Returns a tuple of (WT SASA, MUT SASA)"""
    # get all the file names that end with .asa.pdb
    file_names = [file for file in os.listdir(variant_path) if file.endswith(".asa.pdb")]
    residue_number = tools.get_variant_location_from_variant_folder(tools.get_folder_name_from_path(variant_path))
    # convert the files to dataframes
    if len(file_names) != 2:
        print(f"WARNING: There are {len(file_names)} files that end with .asa.pdb in {variant_path}.")
        return None, None
    for file in file_names:
        if file.startswith('AF'):
            wt_sasa_df = tools.convert_pdb_to_dataframe(os.path.join(variant_path, file))
        else:
            mut_sasa_df = tools.convert_pdb_to_dataframe(os.path.join(variant_path, file))
    # extract the SASA from the dataframes, which is the sum of the 11th column where the residue number =residue_number
    wt_sasa = wt_sasa_df[wt_sasa_df[5] == residue_number][11].sum()
    mut_sasa = mut_sasa_df[mut_sasa_df[5] == residue_number][11].sum()
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
    return wt_stability, mut_stability, mut_stability - wt_stability


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
    # convert file to csv
    af_df = tools.convert_pdb_to_dataframe(af_file_path)

    # get the pLDDT score, which is the first value in the 11th column where the residue number = residue_num
    # in the dataframe
    value = af_df[af_df[5] == residue_num][10].iloc[0]
    return value


def get_average_plddt_score(pdb_file_pwd: str):
    pass


def get_substitution_matrix_value(wt_aa: str, mut_aa: str):
    """Returns the substitution value of the given mutation.
    Uses the BLOSUM62 matrix.
    Args:
        wt_aa: the wild type amino acid
        mut_aa: the mutated amino acid
    Returns:
        int: the substitution value of the given mutation, according to the BLOSUM62 matrix.
    """
    return sm.get_blosum62_value(wt_aa, mut_aa)


def get_consurf_conservation_score_for_all_residues(path_to_gene_folder: str) -> dict:
    """Returns the conservation score of all the residues in the protein.
    Uses the consurf file, which is a text file downloaded from the consurf website.
    The conservation score is the 5th column in the file, and the residue number is in the first column.
    The function finds the value in the 5th column and in the row that corresponds to the given residue number.
    It skips the beginning of the file, which ends with a line that starts with "POS".
    Args:
        path_to_gene_folder: the path to the gene folder, to get the consurf file
    Returns:
        dict: a dictionary of {residue number: conservation score}
    """
    # get the consurf file
    consurf_file = [file for file in os.listdir(path_to_gene_folder) if file.startswith("consurf")][0]
    # read the file
    with open(os.path.join(path_to_gene_folder, consurf_file), 'r') as f:
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


def get_consurf_conservation_score(path_to_gene_folder: str, path_to_variant_folder: str) -> int:
    """Returns the conservation score of the given residue.
    Uses the consurf file, which is a text file downloaded from the consurf website.
    The conservation score is the 5th column in the file, and the residue number is in the first column.
    The function finds the value in the 5th column and in the row that corresponds to the given residue number.
    It skips the beginning of the file, which ends with a line that starts with "POS".
    Args:
        path_to_gene_folder: the path to the gene folder, to get the consurf file
        path_to_variant_folder: the path to the variant folder, to get the residue number
    Returns:
        float: the conservation score of the given residue

    Use:
        consurf = ext_feat.get_consurf_conservation_score("/home/inbar/check/ACTB/", "/home/inbar/check/ACTB/ACTB_P60709_Q189R")
    """
    # get the folder name from the path
    folder_name = tools.get_folder_name_from_path(path_to_variant_folder)
    # get the residue number from the folder name
    residue_num = tools.get_residue_number_from_variant_folder(folder_name)
    residue_num = str(residue_num)  # convert to string
    score = 0  # the score to be returned
    # read the consurf file
    consurf_file_path = os.path.join(path_to_gene_folder, "consurf_grades.txt")
    with open(consurf_file_path, 'r') as f:
        for line in f:
            if line.startswith(" POS"):
                break
        for line in f:
            # find the line that starts with the residue number
            if line.strip().startswith(residue_num):
                # get the conservation score
                score = line.split()[4]
                break

    # if the residue number was not found, write an error message to the log file
    if score == 0:
        with open("log.txt", 'a') as log_file:
            log_file.write(f"Residue number {residue_num} in path {path_to_gene_folder} was not found in the consurf"
                           f" file.\n")
    if "*" in str(score):  # if the score has a star, it means that the conservation score is not reliable
        return 10  # 10 stands for not reliable
    else:
        return int(score)


