import subprocess as sp
import blosum as sm  # blosum matrix, from blosum.py in this directory
import csv
import os
import general_tools as tools
import residue_properties as rp
from biopandas.pdb import PandasPdb


def get_bfactor_from_pdb(pdb_path: str) -> list:
    """Extracts the B-factor from the given PDB file.
    Args:
        pdb_path (str): The path to the PDB file.
    Returns:
        list: A list of the B-factor
    """
    ppdb = PandasPdb()
    ppdb.read_pdb(pdb_path)
    bfactor = ppdb.df['ATOM']['b_factor']
    return bfactor


def extract_feature(feature: str, variant_path: str, wt_or_mut: str = None, aa1: str = None, aa2: str = None) -> float:
    """Extracts the given feature from the given variant.
    Args:
        feature (str): The feature to extract, e.g. 'blosum', 'hydrophobicity', 'volume'.
        variant_path (str): The path to the variant folder.
        wt_or_mut (str): 'wt' or 'mut'.
        aa1 (str): The first amino acid.
        aa2 (str): The second amino acid.
    Returns:
        float: The value of the feature.
    """
    feature_mapping = {
        'blosum': get_substitution_matrix_value,
        'hydrophobicity': get_hydrophobicity,
        'volume': get_aa_volume,
        'plddt': get_plddt,
        'oda': get_oda,
        'opra': get_opra,
        'sasa': get_sasa,
        'stability': get_stability
    }

    is_only_for_wt = ['plddt']
    is_for_wt_and_mut = ['oda', 'opra', 'sasa', 'stability']
    input_both_aa = ['blosum']
    input_aa = ['hydrophobicity', 'volume']

    if feature in feature_mapping:
        feature_function = feature_mapping[feature]
        if feature in input_aa:
            return feature_function(aa1)
        elif feature in input_both_aa:
            return feature_function(aa1, aa2)
        else:
            return feature_function(variant_path)
    else:
        raise ValueError(f"Feature {feature} is not supported.")


def run_oda(variant_path: str, wt_or_mut: str):
    """Creates a file with the ODA score for the given PDB file.
    The file is created in the variant folder.
    Args:
        variant_path (str): The path to the variant folder.
        wt_or_mut (str): 'wt' or 'mut'.
    """
    # change the working directory to the variant folder
    os.chdir(variant_path)

    # get gene_id from folder name
    folder_name = tools.get_folder_name_from_path(variant_path)
    variant_location = tools.get_variant_location_from_variant_folder(folder_name)
    gene_id = folder_name.split('_')[1]

    if wt_or_mut == 'mut':
        amino_acid = tools.get_amino_acid_of_variant_from_variant_folder(folder_name)
        three_letter_amino_acid = tools.convert_1_letter_aa_to_3_letter(amino_acid)
        command = f"pyDock3 {three_letter_amino_acid}{variant_location}_AF_{gene_id}.pdb oda"
    elif wt_or_mut == 'wt':
        command = f"pyDock3 AF_{gene_id}.pdb oda"
    else:
        raise ValueError("wt_or_mut must be 'wt' or 'mut'")
    sp.run(command, shell=True)


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
    # Open the file and get the values
    for file in file_names:
        with open(os.path.join(variant_path, file), 'r') as f:
            lines = f.readlines()
        # get the SASA
        sasa = float(lines[1].split()[1])
        # get the WT and MUT
        if file.startswith("AF"):
            wt_sasa = sasa
        else:
            mut_sasa = sasa
    # get the SASA for the WT and MUT
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


def get_hydrophobicity(aa: str) -> float:
    """Returns the hydrophobicity of the residue"""
    return rp.get_aa_hydrophobicity_by_kd_scale(aa)


def get_aa_volume(aa: str) -> float:
    """Returns the volume of the residue"""
    return rp.get_aa_volume(aa)


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


def get_plddt(residue_num: str, pdb_file_pwd: str):
    """returns the pLDDT score of the residue
    Uses the alphafold pdb file to extract the pLDDT score"""

    # convert file to csv
    pdb_file_pwd = convert_pdb_to_csv(pdb_file_pwd)

    value = None  # The value to be returned
    # read the csv file
    with open(pdb_file_pwd, 'r') as csv_file:
        reader = csv.reader(csv_file)
        for row in reader:
            if row[5] == residue_num:
                value = row[10]
                break  # Exit the loop after finding the first row
    if value is None:
        raise ValueError("The residue number is not valid, or the residue is not in the pdb file.")

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


