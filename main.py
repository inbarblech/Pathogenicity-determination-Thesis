import shutil

import uniprot_info as uni
import add_mutation as add_mut
import dvd_data as dvd
import os
import extract_features as ext_feat
import logging as log
import pandas as pd
import general_tools as tools
import write_data_to_csv as write_csv

PATH_TO_DATA_FOLDER = "/home/inbar/DVDdata/"
PATH_TO_CSV = "/home/inbar/all_variants.csv"
PATH_TO_VARIANTS_FOLDER = "/home/inbar/variants/"
PATHOGENICITY_TYPES = ["Benign", "Pathogenic"]


def create_variants_folders_from_csv(csv_path, output_folder):
    """Reads the csv and creates folders for each variant."""
    df = pd.read_csv(csv_path)

    # check if pathogenicity folders exist
    for pathogenicity in PATHOGENICITY_TYPES:
        path = f"{output_folder}{pathogenicity}/"
        if not os.path.exists(path):
            create_pathogenicity_folders(output_folder)

    existing_folders = []
    # Create folders for each variant, inside the right pathogenicity folder
    for index, row in df.iterrows():
        gene = row["gene"]
        variant = row["variant"]
        variant1letter = tools.convert_variant_to_1_letter(variant)
        uniprot_id = row["uniprot_id"]
        pathogenicity = row["pathogenicity"]

        # Create the gene folder if it doesn't exist
        gene_folder_path = f'{output_folder}{pathogenicity}/{gene}'
        if not os.path.exists(gene_folder_path):
            os.makedirs(gene_folder_path)
        else:
            print(f"Gene folder already exists: {gene_folder_path}")

        # Create the variant folder inside the gene folder
        variant_folder = f'{gene}_{uniprot_id}_{variant1letter}'
        variant_folder_path = os.path.join(gene_folder_path, variant_folder)
        if not os.path.exists(variant_folder_path):
            os.makedirs(variant_folder_path)
            print(f"Variant folder created: {variant_folder_path} for variant {variant} from gene {gene}")
        else:
            print(f"Variant folder already exists: {variant_folder_path}")
            existing_folders.append(f"variant {variant} from gene {gene}")

    print(f"Existing folders: {existing_folders}")
    return existing_folders


def count_number_of_variant_folders(path_to_folder):
    """Counts the number of folders and sub-folders in the given folder."""
    counter = 0
    for root, dirs, files in os.walk(path_to_folder):
        for dir in dirs:
            counter += 1
    return counter


def create_pathogenicity_folders(path_to_folder):
    """Creates the pathogenicity folders (Benign and Pathogenic) in the given folder."""
    for pathogenicity in PATHOGENICITY_TYPES:
        path = f"{path_to_folder}{pathogenicity}/"
        os.makedirs(path)

# def write_pdb_file_to_folder(variant_folder_path):
#     """Writes the pdb file to the variant folder."""
#     uni_id = variant_folder_path.split("_")[1]
#     # change directory to variant folder
#     os.chdir(variant_folder_path)
#     # get pdb file from alpha fold and write it to the variant folder
#     add_mut.get_structure_af(uni_id)


def write_pdb_files_to_all_variants_in_gene_folder(gene_folder_path: str) -> None:
    """Writes the pdb files to all the variant folders in the gene folder.
    Uses the uniprot id to get the pdb file from alpha fold.
    Args:
        gene_folder_path: path to the gene folder.
    """
    # change directory to gene folder
    os.chdir(gene_folder_path)
    gene_name = gene_folder_path.split("/")[-1] # get gene name from path, e.g. "BRCA1" from "path/to/BRCA1"
    uni_id = uni.get_uniprot_id(gene_name)
    add_mut.get_structure_af(uni_id)


def write_pdb_files_to_all_variants_in_all_gene_folders(path_to_gene_folders: str) -> None:
    """Writes the pdb files to all the variant folders in all the gene folders.
    Calls write_pdb_files_to_all_variants_in_gene_folder for each gene folder.
    Args:
        path_to_gene_folders: path to the folder containing all the gene folders.
    """
    # change directory to the folder containing all the gene folders
    os.chdir(path_to_gene_folders)
    # get list of gene folders paths
    gene_folders = os.listdir(path_to_gene_folders)
    # get list of gene folder paths
    gene_folder_paths = [os.path.join(path_to_gene_folders, folder) for folder in gene_folders]
    # write pdb files to all variant folders
    for gene_folder in gene_folder_paths:
        write_pdb_files_to_all_variants_in_gene_folder(gene_folder)


def copy_pdb_files_to_all_variant_folders(path_to_gene_folder: str, uniprot_id: str, gene_name) -> None:
    """Creates a copy of the pdb file in the gene folder, in all the variant folders of the gene."""
    errors = []
    # change directory to the folder containing all the gene folders
    os.chdir(path_to_gene_folder)
    # get list of gene folders paths
    variant_folders = os.listdir(path_to_gene_folder)
    # filter out non-folder entries
    variant_folders = [
        folder
        for folder in variant_folders
        if os.path.isdir(os.path.join(path_to_gene_folder, folder)) and len(os.listdir(os.path.join(path_to_gene_folder,
                                                                                                    folder))) == 0]
    # get list of gene folder paths
    gene_folder_paths = [os.path.join(path_to_gene_folder, folder) for folder in variant_folders]
    # copy pdb files to all variant folders
    for variant_folder in gene_folder_paths:
        print(f"Copying pdb file to variant folder {variant_folder}")
        variant = variant_folder.split("_")[2]  # get variant from path, e.g. "A123B" from "path/to/BRCA1_A123B"
        variant_aa = variant[0]  # get the amino acid from the variant, e.g. "A" from "A123B"
        pos = int(variant[1:-1])  # get the position from the variant, e.g. 123 from "A123B"
        seq = uni.get_sequence(gene_name)
        if uni.get_sequence_length(gene_name) < int(pos):  # check if the position is in the sequence
            print(f"Variant {variant} not in sequence for uniprot id {uniprot_id}")
            errors.append(f"Variant {variant} not in sequence for uniprot id {uniprot_id}")
            continue
        if add_mut.check_aa_in_sequence(variant_aa, seq, pos) is False:
            print(f"Variant {variant} not in sequence for uniprot id {uniprot_id}")
            errors.append(f"Variant {variant} not in sequence for uniprot id {uniprot_id}")
            continue
        # Copy the pdb file to the variant folder
        if os.path.isfile(f"{path_to_gene_folder}/AF-{uniprot_id}-F1-model_v4.pdb"):
            shutil.copyfile(f"{path_to_gene_folder}/AF-{uniprot_id}-F1-model_v4.pdb", f"{variant_folder}/AF_{uniprot_id}.pdb")
    # write errors to file in gene folder
    with open(f"{path_to_gene_folder}/errors.txt", "w") as f:
        for error in errors:
            f.write(error + "\n")


def copy_pdb_files_to_all_variant_folders_in_all_gene_folders(path_to_gene_folders: str) -> None:
    """Creates a copy of the pdb file in the gene folder, in all the variant folders of all the genes.
    Calls copy_pdb_files_to_all_variant_folders for each gene folder.
    Args:
        path_to_gene_folders: path to the folder containing all the gene folders.
    """
    # change directory to the folder containing all the gene folders
    os.chdir(path_to_gene_folders)
    # get list of gene folders paths
    gene_folders = os.listdir(path_to_gene_folders)
    # get list of gene folder paths
    gene_folder_paths = [os.path.join(path_to_gene_folders, folder) for folder in gene_folders]
    # copy pdb files to all variant folders
    for gene_folder in gene_folder_paths:
        # get uniprot id from gene folder name
        gene_name = gene_folder.split("/")[-1]  # get gene name from path, e.g. "BRCA1" from "path/to/BRCA1"
        uniprot_id = uni.get_uniprot_id(gene_name)
        copy_pdb_files_to_all_variant_folders(gene_folder, uniprot_id, gene_name)


def extract_features(variant, uni_id, residue_num, pdb_file_path):
    """Extract features"""

    # Extract plddt value
    plddt_residue_value = ext_feat.get_plddt(residue_num, pdb_file_path)  # Get plddt value for residue

    # ODA
    oda_wt = ext_feat.get_oda(variant, uni_id, "wt")
    oda_mut = ext_feat.get_oda(variant, uni_id, "mut")
    delta_oda = ext_feat.get_delta_oda(variant, uni_id)

    # OPRA
    opra_wt = ext_feat.get_opra(variant, uni_id, "wt")
    opra_mut = ext_feat.get_opra(variant, uni_id, "mut")
    delta_opra = ext_feat.get_delta_opra(variant, uni_id)

    # Hydrogen bonds
    hbonds_wt = ext_feat.get_hbonds(variant, uni_id, "wt")
    hbonds_mut = ext_feat.get_hbonds(variant, uni_id, "mut")
    delta_hbonds = ext_feat.get_delta_hbonds(variant, uni_id)

    aa1, aa2 = tools.get_aa1_aa2(variant)

    # Substitution matrix value
    sub_value = ext_feat.get_sub_mat(aa1, aa2)

    # Save features as csv row
    features = {"Gene_name": gene_name, "Variant": variant, "UniProt_ID": uni_id,
                "PLDDT": plddt_residue_value, "ODA_wt": oda_wt, "ODA_mut": oda_mut,
                "Delta_ODA": delta_oda, "OPRA_wt": opra_wt, "OPRA_mut": opra_mut,
                "Delta_OPRA": delta_opra, "HBonds_wt": hbonds_wt, "HBonds_mut": hbonds_mut,
                "Delta_HBonds": delta_hbonds, "Substitution_matrix": sub_value}
    df_variant_row = save_data_as_df(pathogenicity, features)
    df = df.append(df_variant_row)


def main_automation_set_up():  #TODO: Create the functions for this
    """This function is used to set up the automation process for the first time in a new environment.
    It creates the folder structure and downloads the necessary files."""
    # Create folder structure
    create_folder_structure(PATH_TO_VARIANTS_FOLDER)
    # Download files
    download_files(PATH_TO_VARIANTS_FOLDER)


if __name__ == '__main__':
    copy_pdb_files_to_all_variant_folders_in_all_gene_folders(f"{PATH_TO_VARIANTS_FOLDER}Pathogenic/")
    copy_pdb_files_to_all_variant_folders_in_all_gene_folders(f"{PATH_TO_VARIANTS_FOLDER}Benign/")
