import subprocess as sp
import blosum as sm  # blosum matrix, from blosum.py in this directory
import csv
import os


def get_oda(gene_id: str, wt_or_mut: str, amino_acid: str):
    """Creates a file with the ODA score for the given PDB file."""
    if wt_or_mut == 'wt':
        command = f"pyDock3 AF-{gene_id}-F1-model_v4.pdb oda"
    elif wt_or_mut == 'mut':
        command = f"pyDock3 {amino_acid}_AF-{gene_id}-F1-model_v4.pdb oda"
    else:
        raise ValueError("wt_or_mut must be 'wt' or 'mut'")
    sp.run(command, shell=True)


def get_opra(gene_id: str, wt_or_mut: str, amino_acid: str):
    """creates a file with the OPRA score for the given PDB file."""
    if wt_or_mut == 'mut':
        command = f"pyDock3 {amino_acid}_AF-{gene_id}-F1-model_v4.pdb opra"
    elif wt_or_mut == 'wt':
        command = f"pyDock3 AF-{gene_id}-F1-model_v4.pdb opra"
    else:
        raise ValueError("wt_or_mut must be 'wt' or 'mut'")
    sp.run(command, shell=True)


def get_sasa(file_name: str, wt_or_mut: str):
    """Creates a file with the SASA score for the given PDB file."""
    command = f"dr_sasa -m 0 -i {file_name}.{wt_or_mut}.pdb"
    sp.run(command, shell=True)


def get_stability(folder_name: str):
    """Returns the stability of the residue using foldX energetics in the file named 'energies.log', inside the folder
    of the variant.
    The stability is in the"""


def get_hydrophobicity():
    """Returns the hydrophobicity of the residue"""


def get_hydrogen_bonds():
    """Returns the number of hydrogen bonds of residue"""


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


def get_substitution_value(wt_aa: str, mut_aa: str):
    """Returns the substitution value of the given mutation.
    Uses the BLOSUM62 matrix.
    Args:
        wt_aa: the wild type amino acid
        mut_aa: the mutated amino acid
    Returns:
        int: the substitution value of the given mutation, according to the BLOSUM62 matrix.
    """
    return sm.get_blosum62_value(wt_aa, mut_aa)


if __name__ == "__main__":
    # Call get_plddt with the path to the pdb file and the residue number
    print(get_plddt("9", '/home/inbar/variants/Benign/POLR1B/POLR1B_Q9H9Y6_P3S/AF-Q9H9Y6-F1-model_v4.pdb'))
    # Call get_substitution_value with the wild type amino acid and the mutated amino acid
    print(get_substitution_value("A", "C"))
    print(convert_pdb_to_csv('/home/inbar/variants/Benign/POLR1B/POLR1B_Q9H9Y6_P3S/AF-Q9H9Y6-F1-model_v4.pdb'))

