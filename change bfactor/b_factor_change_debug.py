import general_tools as gt
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from data_retrievel_and_feature_extraction import uniprot_info as uni
from Bio.PDB import PDBIO, Model, Chain, Residue, Atom
import pandas as pd
from biopandas.pdb import PandasPdb
from Bio.PDB import PDBParser


def pdb_to_dataframe(pdb_file):
    # Create a PDB parser and load the structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    # Initialize lists to store PDB data
    atom_data = []

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom_record = {
                        'ATOM': atom.get_serial_number(),
                        'NAME': atom.get_name(),
                        'RESNAME': residue.get_resname(),
                        'CHAINID': chain.id,
                        'RESIDUESEQ': residue.id[1],
                        'X': atom.get_coord()[0],
                        'Y': atom.get_coord()[1],
                        'Z': atom.get_coord()[2],
                        'OCCUPANCY': atom.get_occupancy(),
                        'TEMPFACTOR': atom.get_bfactor(),
                        'ELEMENT': atom.element,
                    }
                    atom_data.append(atom_record)

    # Create a DataFrame from the atom data
    df = pd.DataFrame(atom_data)

    return df


######### Functions for addition of b_factor column by pathogenicity to pdb_df dataframe #########


def create_dictionary_with_all_residues(protein_length):
    """
    Creates a dictionary where the keys are residue numbers and the values are
    the residue names. The dictionary contains all residues from 1 to the
    protein_length.
    The values for all residues are set to 'None' by default.
    """
    # Create an empty dictionary
    residue_dict = {}
    # Iterate through the range of numbers from 1 to the protein_length
    for residue_number in range(1, protein_length + 1):
        # Add the residue number and residue name to the residue_dict dictionary
        residue_dict[residue_number] = None
    return residue_dict


def get_b_factor_dict(residue_dict, pathogenicity_df, gene_length):
    """
    Go over all the residues: Create a subset of the pathogenicity_df by residue.
    For each subset, if there is only one value, then insert the pathogenicity value into the dictionary.
    If there are multiple values, then check if they are the same. If they are, then insert the pathogenicity value into the dictionary.
    If they are not, then insert 'Vague' into the dictionary.
    """
    for residue in range(gene_length):
        # Create a subset of the pathogenicity_df by residue
        pathogenicity_by_residue = pathogenicity_df[pathogenicity_df['position'] == residue]
        # Check if there is only one value in the pathogenicity_by_residue dataframe
        if len(pathogenicity_by_residue) == 1:
            # Get the pathogenicity value from the pathogenicity_by_residue dataframe
            pathogenicity = pathogenicity_by_residue['pathogenicity'].values[0]
            # Insert the pathogenicity value into the residue_dict dictionary
            residue_dict[residue] = pathogenicity
        elif len(pathogenicity_by_residue) > 1:
            # Get the pathogenicity values from the pathogenicity_by_residue dataframe
            pathogenicity_values = pathogenicity_by_residue['pathogenicity'].values
            # Check if all the pathogenicity values are the same
            if all(pathogenicity == pathogenicity_values[0] for pathogenicity in pathogenicity_values):
                # Insert the pathogenicity value into the residue_dict dictionary
                residue_dict[residue] = pathogenicity_values[0]
            else:
                # Insert 2 for Vague information into the residue_dict dictionary
                residue_dict[residue] = 2
    return residue_dict


def add_b_factor_by_residue(pdb_df, b_factor_dict):
    """
    Adds a b_factor column to the pdb_df dataframe. The b_factor column
    contains the b_factor values for each residue in the pdb_df dataframe.
    The b_factor values are taken from the b_factor_dict dictionary.
    Addition is done by residue, so all atoms in a residue will have the same b_factor.
    """
    # Make sure there's a bfactor column in the pdb_df dataframe
    if 'TEMPFACTOR' not in pdb_df.columns:
        # raise an error if there's no bfactor column
        raise ValueError('There is no bfactor column in the pdb_df dataframe')
    # Iterate through the rows of the pdb_df dataframe
    for index, row in pdb_df.iterrows():
        # Get the residue number for the current row
        residue_number = row['RESIDUESEQ']
        # Get the b_factor value for the current row
        b_factor = b_factor_dict[residue_number]
        # Set the b_factor value for all rows with the current residue number
        pdb_df.loc[pdb_df['RESIDUESEQ'] == residue_number, 'TEMPFACTOR'] = b_factor
    return pdb_df


def dataframe_to_pdb(dataframe, output_file):
    ppdb = PandasPdb()
    ppdb.df['ATOM'] = dataframe
    ppdb.to_pdb(output_file, records=['ATOM'])
#%%


# Define a function to update the B-factor (11th column) value
def update_bfactor(line, new_value):
    line_list = list(line)
    line_list[60:70] = f'{new_value:6.2f}'.encode('utf-8')
    return bytes(line_list)


def create_pdb_with_bfactor(pdb_file):
    # Path to your input and output files
    input_file_path = pdb_file
    output_file_path = "COL4A5_with_b_factor.pdb"

    # Flag to indicate whether you're inside the "ATOM" table
    inside_atom_table = False

    # Read the input file and write to the output file
    with open(input_file_path, "r") as input_file, open(output_file_path, "wb") as output_file:
        for line in input_file:
            # Check if the line starts with "ATOM" to identify the start of the table
            if line.startswith("ATOM"):
                inside_atom_table = True

            if inside_atom_table:
                # Perform your calculations here
                # For example, let's increment the B-factor by 10 for demonstration
                bfactor = float(line[60:70])
                new_bfactor = bfactor + 10.0
                line = update_bfactor(line, new_bfactor)

            # Write the line to the output file
            output_file.write(line.encode('utf-8'))


if __name__ == '__main__':

    genes = ["COL2A1", "COL4A3", "COL4A5", "WFS1", "SLC26A4", "MYO7A", "FGFR1", "GJB2"]
    gene = genes[7]
    pdb_file = f'C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\AF_structures\\{gene}.pdb'

    pdb_df = pdb_to_dataframe(pdb_file)
    print(pdb_df)
    print(pdb_df.columns)

    gene_length = uni.get_sequence_length(gene)
    pathogenicity_df = pd.read_csv(f'C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\gene_specific_df\\{gene}_combined_with_source_and_position.csv')
    print(pathogenicity_df)

    # leave only the pathogenicity and residue number columns
    pathogenicity_df = pathogenicity_df[['position', 'pathogenicity']]
    # Change pathogenicity to 1 and 0, where 1 is pathogenic and 0 is benign
    pathogenicity_df['pathogenicity'] = pathogenicity_df['pathogenicity'].map({'pathogenic': 1, 'benign': 0})

    residue_dict = create_dictionary_with_all_residues(gene_length)

    b_factor_dict = get_b_factor_dict(residue_dict, pathogenicity_df, gene_length)
    pdb_df = add_b_factor_by_residue(pdb_df, b_factor_dict)
    print(b_factor_dict)

    benign_residues = 0
    pathogenic_residues = 0
    vague_residues = 0
    no_data_residues = 0

    for key in b_factor_dict:
        if b_factor_dict[key] == 0:
            benign_residues += 1
        elif b_factor_dict[key] == 1:
            pathogenic_residues += 1
        elif b_factor_dict[key] == 2:
            vague_residues += 1
        elif b_factor_dict[key] == None:
            no_data_residues += 1

    print(f"Data for gene: {gene}")
    print(f'benign_residues: {benign_residues}')
    print(f'pathogenic_residues: {pathogenic_residues}')
    print(f'vague_residues: {vague_residues}')
    print(f'no_data_residues: {no_data_residues}')

    percentage_benign = benign_residues / gene_length
    percentage_pathogenic = pathogenic_residues / gene_length
    percentage_vague = vague_residues / gene_length
    percentage_no_data = no_data_residues / gene_length

    # Create a graph showing the percentage of each value in the bfactor dictionary
    labels = ['Benign', 'Pathogenic', 'Mixed', 'No data']
    sizes = [percentage_benign, percentage_pathogenic, percentage_vague, percentage_no_data]
    colors = ['gold', 'yellowgreen', 'lightcoral', 'lightskyblue']
    explode = (0.1, 0, 0, 0)  # explode 1st slice

    plt.pie(sizes, labels=labels, colors=colors,
            autopct='%1.1f%%', shadow=True, startangle=140)

    plt.title(f'Percentage of residues of type in model data for {gene}')
    # add title for each value
    # plt.legend(labels, loc="best")
    plt.axis('equal')
    plt.show()
