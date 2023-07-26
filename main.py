import csv
import shutil

import subprocess as sp
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
PATH_TO_CSV_BENIGN = "/home/inbar/all_BENIGN_variants.csv"
PATH_TO_CSV_PATHOGENIC = "/home/inbar/all_PATHOGENIC_variants.csv"
PATH_TO_VARIANTS_FOLDER = "/home/inbar/variants/"
PATHOGENICITY_TYPES = ["Benign", "Pathogenic"]


def create_variants_folders_from_csv(csv_path, output_folder, pathogenicity):
    """Reads the csv and creates folders for each variant.
    Args:
        csv_path (str): The path to the csv file.
        output_folder (str): The path to the output folder.
        pathogenicity (str): The pathogenicity type, either "Benign" or "Pathogenic".
        """
    df = pd.read_csv(csv_path)

    # # check if pathogenicity folders exist
    # for pathogenicity in PATHOGENICITY_TYPES:
    #     path = f"{output_folder}{pathogenicity}/"
    #     if not os.path.exists(path):
    #         create_pathogenicity_folders(output_folder)

    existing_folders = []
    # Create folders for each variant, inside the right pathogenicity folder
    for index, row in df.iterrows():
        gene = row["gene"]
        variant = row["variant"]
        variant1letter = tools.convert_variant_to_1_letter(variant)
        uniprot_id = row["uniprot_id"]

        # Create the gene folder if it doesn't exist
        gene_folder_path = f'{output_folder}/{gene}'
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


def extract_stability_data_to_csv(path_to_pathogenicity_folder: str):
    """Uses ext_feat.get_stability to extract the stability data from the files in the genes folders."""
    # change directory to the pathogenicity folder
    os.chdir(path_to_pathogenicity_folder)
    # get list of gene folders paths
    gene_folders = os.listdir(path_to_pathogenicity_folder)
    # get list of gene folder paths
    gene_folder_paths = [os.path.join(path_to_pathogenicity_folder, folder) for folder in gene_folders]
    # extract stability data from all variant folders
    with open("stability_data.csv", "w") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["gene", "variant", "wt", "mut", "delta"])

        for gene_folder in gene_folder_paths:
            if not os.path.isdir(gene_folder):
                continue
            variant_folders = os.listdir(gene_folder)
            variant_folder_paths = [os.path.join(gene_folder, folder) for folder in variant_folders]
            gene = gene_folder.split("/")[-1]
            for variant_folder in variant_folder_paths:
                if not os.path.isdir(variant_folder):
                    continue
                # Check if exists a file that ends with "scanning_output.txt"
                if any([file.startswith("energies") for file in os.listdir(variant_folder)]):
                    variant = variant_folder.split("/")[-1].split("_")[-1]
                    wt, mut, delta = ext_feat.get_stability(variant_folder)
                    writer.writerow([gene, variant, wt, mut, delta])
                    print(f"Stability data extracted for variant {variant} from gene {gene}, it is {delta}")


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
    variant_folder_paths = [os.path.join(path_to_gene_folder, folder) for folder in variant_folders]
    # copy pdb files to all variant folders
    for variant_folder in variant_folder_paths:
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


def create_foldx_mutation_to_variants_in_gene_folder(path_to_gene_folder: str) -> None:
    """This function creates a FoldX mutation file for all the variants in a gene folder.
    Args:
        path_to_gene_folder: path to the gene folder, that contains all the variant folders.
        uniprot_id: uniprot id of the gene.
        gene_name: name of the gene.
    """
    print("entered create_foldx_mutation_to_variants_in_gene_folder function")
    errors = []
    # change directory to the folder containing all the gene folders
    os.chdir(path_to_gene_folder)
    # get list of variant folders paths
    all_items_in_gene_folder = os.listdir(path_to_gene_folder)

    # get list of variant folder paths
    variant_folders_paths = [os.path.join(path_to_gene_folder, item) for item in all_items_in_gene_folder if
                            os.path.isdir(os.path.join(path_to_gene_folder, item))]

    print(f"variant_folders_paths: {variant_folders_paths}")
    # create mutation to all variants in all genes, using the function add_mut.create_mut_using_foldx
    for variant_folder in variant_folders_paths:
        folder_name = os.path.basename(variant_folder)
        uniprot_id = folder_name.split("_")[1]
        mut = folder_name.split("_")[2]
        if not os.path.isfile(f"{variant_folder}/AF_{uniprot_id}.pdb"):
            print(f"No alphafold pdb file in variant folder {variant_folder}")
            errors.append(f"No alphafold pdb file in variant folder {variant_folder}")
            continue
        os.chdir(variant_folder)
        add_mut.create_mutation_using_foldx(mut, uniprot_id)
    # write errors to file in gene folder
    with open(f"{path_to_gene_folder}/errors.txt", "w") as f:
        for error in errors:
            f.write(error + "\n")


def create_mutation_file_to_all_variants(path: str) -> None:
    """Uses the function add_mut.create_mutation to create a mutation to all variants in all genes in the folder.
    Args:
        path: path to the folder containing all the gene folders.
    """
    errors = []
    # change directory to the folder containing all the gene folders
    os.chdir(path)
    # get list of gene folders paths
    gene_folders = os.listdir(path)
    # get list of gene folder paths, do not include files
    gene_folder_paths = [os.path.join(path, folder) for folder in gene_folders if os.path.isdir(os.path.join(path, folder))]
    # create mutation to all variants in all genes
    for gene_folder in gene_folder_paths:
        print(gene_folder)
        if not os.path.isfile(f"{gene_folder}/errors.txt"):
            print(f"No errors.txt file in gene folder {gene_folder}")
            errors.append(f"No errors.txt file in gene folder {gene_folder}")
            continue
        # Check if there is a file of type pdb in the gene folder
        if not any(file.endswith(".pdb") for file in os.listdir(gene_folder) if
                   os.path.isfile(os.path.join(gene_folder, file))):
            print(f"No pdb file in gene folder {gene_folder}")
            errors.append(f"No pdb file in gene folder {gene_folder}")
            continue

        create_foldx_mutation_to_variants_in_gene_folder(f"{gene_folder}/")
        # return to the path containing all the gene folders


def extract_features(path: str, output_file_path: str) -> list:
    """Extract features for all variants in all genes, and write to csv file.
    The function checks if there is a pdb file in the variant folder. If there is, it extracts the features.
    If there is not, it skips the variant and adds the variant and the path to the errors file.
    Args:
        path: path to the folder containing all the gene folders.
        output_file_path: path to the output file (csv file).
    returns:
        errors: list of errors
    """
    os.chdir(path)
    # get list of gene folders paths
    gene_folders = os.listdir(path)
    # get list of gene folder paths
    gene_folder_paths = [os.path.join(path, folder) for folder in gene_folders]
    # create list of features
    features = []
    # create list of errors
    errors = []
    # extract features for all variants in all genes
    for gene_folder in gene_folder_paths:
        # get uniprot id from gene folder name
        gene_name = gene_folder.split("/")[-1]  # get gene name from path, e.g. "BRCA1" from "path/to/BRCA1"
        uni_id = uni.get_uniprot_id(gene_name)
        # get variant from gene folder name, e.g. "A123B" from "KCNQ1_P51787_A2V" then variant="A2V"
        mut = gene_folder.split("_")[-1]
        residue_num = mut[1:-1]
        # get path to pdb file
        pdb_file_path = f"{gene_folder}/AF_{uni_id}.pdb"
        extract_features_per_variant(mut, uni_id, residue_num, pdb_file_path)
    return errors


def extract_features_per_variant(variant, uni_id, residue_num, pdb_file_path) -> list:
    """Extract features
    Args:
        variant: variant, e.g. "A2V"
        uni_id: uniprot id, e.g. "P51787"
        residue_num: residue number, e.g. "123"
        pdb_file_path: path to pdb file, e.g. "path/to/KCNQ1_P51787_A2V/AF_P51787.pdb"
    returns:
        errors: list of errors
        """

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
    sub_value = ext_feat.get_substitution_matrix_value(aa1, aa2)

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

    # create foldX mutation file to all variants in all gene folders
    create_mutation_file_to_all_variants(f"{PATH_TO_VARIANTS_FOLDER}Pathogenic/")
    create_mutation_file_to_all_variants(f"{PATH_TO_VARIANTS_FOLDER}Benign/")

#
# if __name__ == '__main__':
#     copy_pdb_files_to_all_variant_folders_in_all_gene_folders(f"{PATH_TO_VARIANTS_FOLDER}Pathogenic/")
#     copy_pdb_files_to_all_variant_folders_in_all_gene_folders(f"{PATH_TO_VARIANTS_FOLDER}Benign/")

if __name__ == "__main__":
    # (wt, mut, delta) = ext_feat.get_stability("/home/inbar/check/ANKH/ANKH_Q9HCJ1_R453W/")
    # print(wt, mut, delta)
    # consurf = ext_feat.get_consurf_conservation_score("/home/inbar/check/ACTB/", "/home/inbar/check/ACTB/ACTB_P60709_Q189R")
    # print(consurf)
    # # ext_feat.get_consurf_conservation_score_for_all_residues("/home/inbar/check/ACTB/")
    # ext_feat.get_sasa("/home/inbar/check/ACTB/ACTB_P60709_V209M", "wt")
    # ext_feat.get_sasa("/home/inbar/check/ACTB/ACTB_P60709_V209M", "mut")
    sp.run("pyDock3", shell=True)