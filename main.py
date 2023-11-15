import csv
import shutil

from data_retrievel_and_feature_extraction import uniprot_info as uni
from data_retrievel_and_feature_extraction import add_mutation as add_mut
import os
from data_retrievel_and_feature_extraction import feature_extraction_column_by_column as ext_feat
import pandas as pd
import general_tools as tools
from data_retrievel_and_feature_extraction import run_command_line_programs_to_folders as run_command_line_programs
import matplotlib.pyplot as plt
from data_retrievel_and_feature_extraction import write_to_folders as write_to_folders

PATH_TO_DATA_FOLDER = "/home/inbar/DVDdata/"
PATH_TO_CSV_BENIGN = "/home/inbar/all_BENIGN_variants.csv"
PATH_TO_CSV_PATHOGENIC = "/home/inbar/all_PATHOGENIC_variants.csv"
PATH_TO_VARIANTS_FOLDER = "/home/inbar/variants/"
PATH_TO_OUTPUT_FOLDER = "/home/inbar/results/"
PATH_TO_INBAR = "/home/inbar/"
PATHOGENICITY_TYPES = ["Benign", "Pathogenic"]
EXAMPLE_VARIANT_PATH = "/home/inbar/variants/Benign/AIFM1/AIFM1_O95831_E92K"
PATH_TO_FEATURES_CSV = "/home/inbar/results/features.csv"
PATH_TO_SMALL_FEATURES_CSV = "/home/inbar/results/small_features.csv"
PATH_TO_NOVEL_BENIGN_VARIANTS = "/home/inbar/variants"


def run_on_gene_folders(path, list_of_genes_to_run_on=None):
    """Runs a function on all gene folders in path.
    Usage:
    run_on_gene_folders(PATH_TO_VARIANTS_FOLDER)"""
    errors = []
    os.chdir(path)
    # get list of gene folders paths
    if list_of_genes_to_run_on is None:
        items_in_folder = os.listdir(path)
    else:
        items_in_folder = list_of_genes_to_run_on
    gene_folders_paths = [os.path.join(path, item) for item in items_in_folder if
                          os.path.isdir(os.path.join(path, item))]
    for gene_folder in gene_folders_paths:
        error = run_command_line_programs_on_folder(gene_folder)
        errors.append(error)

    with open(f"{path}/log_oda_opra_sasa.txt", "w") as f:
        for error in errors:
            f.write(f"{error}\n")


def run_command_line_programs_on_folder(path_to_gene_folder):
    """Runs command line programs on all variant folders in path_to_folder"""
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
    already_ran = True
    for variant_folder in variant_folders_paths:
        # check if there is already a file of type oda in the folder
        if already_ran:
            oda_files = [file for file in os.listdir(variant_folder) if file.endswith(".oda")]
            sasa_files = [file for file in os.listdir(variant_folder) if file.endswith(".asa")]
            if len(oda_files) == 0 or len(sasa_files) == 0:
                print("Starting to run oda and sasa")
                already_ran = False
            else:
                continue
        try:
            run_command_line_programs.run_opra_oda_sasa(variant_folder)
        except Exception as e:
            errors.append(e)

    with open(f"{path_to_gene_folder}/log_oda_opra_sasa_mut.txt", "w") as f:
        for error in errors:
            f.write(f"{error}\n")
    return errors


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


def copy_pdb_files_to_all_variant_folders(path_to_gene_folder: str, uniprot_id: str = None, gene_name: str = None) -> None:
    """Creates a copy of the pdb file in the gene folder, in all the variant folders of the gene."""
    if gene_name is None:
        gene_name = path_to_gene_folder.split("/")[-1]  # get gene name from path, e.g. "BRCA1" from "path/to/BRCA1"
        print(gene_name)
    if uniprot_id is None:
        uniprot_id = uni.get_uniprot_id(gene_name)
    print(uniprot_id)
    print(gene_name)
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
        variant_folder_name = variant_folder.split("/")[-1]  # get variant folder name from path, e.g. "BRCA1_A123B" from "path/to/BRCA1_A123B"
        variant = variant_folder_name.split("_")[2]  # get variant from path, e.g. "A123B" from "path/to/BRCA1_A123B"
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
        if os.path.isfile(f"{path_to_gene_folder}/AF_{uniprot_id}.pdb"):
            shutil.copyfile(f"{path_to_gene_folder}/AF_{uniprot_id}.pdb", f"{variant_folder}/AF_{uniprot_id}.pdb")
            print(f"File {path_to_gene_folder}/AF_{uniprot_id}.pdb copied to {variant_folder}/AF_{uniprot_id}.pdb")
        else:
            print(f"File {path_to_gene_folder}/AF_{uniprot_id}.pdb does not exist")
            errors.append(f"File {path_to_gene_folder}/AF_{uniprot_id}.pdb does not exist")
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


def create_foldx_mutation_to_variants_in_gene_folder(path_to_gene_folder: str, count=None) -> None:
    """This function creates a FoldX mutation file for all the variants in a gene folder.
    Args:
        path_to_gene_folder: path to the gene folder, that contains all the variant folders.
        uniprot_id: uniprot id of the gene.
        gene_name: name of the gene.
    """
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
        foldx_mutation_exists = False
        files = os.listdir(variant_folder)
        for file in files:
            if file.endswith(".txt"):
                print(f"Mutation file already exists for variant folder {variant_folder}")
                foldx_mutation_exists = True
                break
        if foldx_mutation_exists:
            continue
        if count is not None:
            count += 1
        print(f"Creating mutation file for variant folder {variant_folder}")
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
    Usage:
        create_mutation_file_to_all_variants(f"{PATH_TO_VARIANTS_FOLDER}Benign_for_gene_specific/
    """
    errors = []
    # change directory to the folder containing all the gene folders
    os.chdir(path)
    # get list of gene folders paths
    gene_folders = os.listdir(path)
    # get list of gene folder paths, do not include files
    gene_folder_paths = [os.path.join(path, folder) for folder in gene_folders if os.path.isdir(os.path.join(path, folder))]
    # create mutation to all variants in all genes
    count = 0
    for gene_folder in gene_folder_paths:
        print(gene_folder)
        # if not os.path.isfile(f"{gene_folder}/errors.txt"):
        #     print(f"No errors.txt file in gene folder {gene_folder}")
        #     errors.append(f"No errors.txt file in gene folder {gene_folder}")
        #     continue
        # Check if there is a file of type pdb in the gene folder
        if not any(file.endswith(".pdb") for file in os.listdir(gene_folder) if
                   os.path.isfile(os.path.join(gene_folder, file))):
            print(f"No pdb file in gene folder {gene_folder}")
            errors.append(f"No pdb file in gene folder {gene_folder}")
            continue

        create_foldx_mutation_to_variants_in_gene_folder(f"{gene_folder}/", count)
        # return to the path containing all the gene folders


def main_automation_set_up():  #TODO: Create the functions for this
    """This function is used to set up the automation process for the first time in a new environment.
    It creates the folder structure and downloads the necessary files."""
    # Create folder structure
    create_folder_structure(PATH_TO_VARIANTS_FOLDER)
    # Download files
    download_files(PATH_TO_VARIANTS_FOLDER)

    # create foldX mutation file to all variants in all gene folders
    create_mutation_file_to_all_variants(f"{PATH_TO_VARIANTS_FOLDER}Benign_for_gene_specific/")
    # create_mutation_file_to_all_variants(f"{PATH_TO_VARIANTS_FOLDER}Benign/")


def get_dataframe_of_number_of_variants_per_gene_per_pathogenicity(features_df) -> pd.DataFrame:
    """Creates a dataframe that contains the number of variants in each pathogenicity for each gene.
    For example:
    Gene    Pathogenic  Benign
    ADHFK   5           12
    FKLK    1           0
    IOI5    120         9

    Input: Dataframe of all data (features.csv)
    Output: Dataframe of number of variants per pathogenicity per gene.

    Usage:
    df = get_dataframe_of_number_of_variants_per_gene_per_pathogenicity(features_df)
    """

    features_df_path = features_df.loc[features_df['pathogenicity'] == 'pathogenic']
    features_df_benign = features_df.loc[features_df['pathogenicity'] == 'benign']
    num_of_variants_per_gene_dict_path = tools.get_number_of_variants_per_gene_dict(features_df_path)
    num_of_variants_per_gene_dict_benign = tools.get_number_of_variants_per_gene_dict(features_df_benign)
    # create dataframe from the two dictionaries
    # Convert dictionaries to DataFrames
    df_pathogenic_values = pd.DataFrame(list(num_of_variants_per_gene_dict_path.items()),
                                        columns=['gene', 'pathogenic'])
    df_benign_values = pd.DataFrame(list(num_of_variants_per_gene_dict_benign.items()), columns=['gene', 'benign'])
    # Merge DataFrames using an outer join to include all keys from both dictionaries
    merged_df = pd.merge(df_pathogenic_values, df_benign_values, on='gene', how='outer')
    # Fill NaN values with a specific value if needed
    merged_df.fillna(0, inplace=True)  # Replace NaN with 0
    return merged_df


def create_one_dataframe_for_genes(genes_list: list, df_list: list) -> pd.DataFrame:
    """This function creates one dataframe for all the genes in the list.
    It adds to the new dataframe every row from the two dataframes that corresponds with the genes in the list."""
    combined_df = pd.DataFrame()
    for gene in genes_list:
        for dataframe in df_list:
            combined_df = combined_df.append(dataframe.loc[dataframe['gene'] == gene])
    return combined_df


def combine_dataframes(df_list: list) -> pd.DataFrame:
    """This function combines all the dataframes in the list to one dataframe.
    It adds to the new dataframe every row from the two dataframes that corresponds with the genes in the list.
    If the dataframes have the same columns, it will combine them to one dataframe.
    If the dataframes have different columns, it will add the columns to the new dataframe and fill the missing values
    with NaN."""
    if not df_list:
        raise ValueError("The list of dataframes is empty.")
    combined_df = df_list[0]
    for dataframe in df_list[1:]:
        combined_df = combined_df.merge(dataframe, on=None, how="outer")
    return combined_df


def create_sub_dataframe_according_to_gene(df: pd.DataFrame, gene: str) -> pd.DataFrame:
    """This function creates a new dataframe with only a given gene from the original dataframe.
    The gene name is a column in the dataframe.
    Usage:
        seven_genes_df = pd.read_csv(f"{PATH_TO_NOVEL_BENIGN_VARIANTS}/Benign_for_gene_specific/combined.csv", header=0)
        WFS1_df = create_sub_dataframe_according_to_gene(seven_genes_df, "WFS1")
        WFS1_df.to_csv("/home/inbar/results/gene_specific_df/WFS1.csv", index=False)
    """
    return df.loc[df['gene'] == gene]


def add_source_column(combined_df, features_df, benign_df) -> pd.DataFrame:
    """This function adds a column to the "combined_df" dataframe that indicates the source of the variant.
    The source is either DVD (for variants from the features_df) or patmut (for variants from the benign_df).
    The check is done by the columns "variant" and "gene".
    """
    source_dict = {}

    for index, row in combined_df.iterrows():
        variant = row['variant']
        gene = row['gene']

        if (variant, gene) in zip(features_df['variant'], features_df['gene']):
            source_dict[index] = "DVD"
        elif (variant, gene) in zip(benign_df['variant'], benign_df['gene']):
            source_dict[index] = "patmut"

    source_df = pd.DataFrame.from_dict(source_dict, orient='index', columns=['source'])
    result_df = combined_df.join(source_df, how='left')

    return result_df


def add_position_column(df: pd.DataFrame) -> pd.DataFrame:
    """This function recieves a dataframe with a column called "variant" and adds a column called "position",
    based on the variant"""
    position_dict = {}
    for index, row in df.iterrows():
        variant = row['variant']
        # get position, which is the number after the first letter
        position = variant[1:]
        # remove the last letter, which is the variant
        position = position[:-1]
        position_dict[index] = position
    position_df = pd.DataFrame.from_dict(position_dict, orient='index', columns=['position'])
    result_df = df.join(position_df, how='left')
    return result_df


def create_predictions_vs_real_csv(path):
    """This function creates a csv file with the predictions and the real values.
    It combines all the csv files (columns position, predictions and pathogenicity) in the path folder to one csv file.
    The file is saved in the path folder.
    Args:
        path: the path to the folder with the predictions csv files.

    Usage:
        path = ("predictions_vs_real/MYO7A")
        create_predictions_vs_real_csv(path)
    """
    # get all the csv files in the path folder
    csv_files = [f for f in os.listdir(path) if f.endswith('.csv') and not f.startswith('predictions_vs_real')]
    # create a dataframe for each csv file
    df_list = []
    for csv_file in csv_files:
        df_list.append(pd.read_csv(f"{path}/{csv_file}", header=0))
    # combine all the dataframes to one dataframe
    combined_df = combine_dataframes(df_list)
    # Remove all columns exceps position, variant, predictions and pathogenicity
    combined_df = combined_df[['position', 'pathogenicity', 'predictions', 'variant']]
    # save the dataframe to a csv file
    combined_df.to_csv(f"{path}/predictions_vs_real_with_variant.csv", index=False)


def count_position_repetitions(df: pd.DataFrame) -> dict:
    """
    This function returns a dictionary where the keys are number of position repetitions and the values are a count of
    how many times this number of repetitions appears in the dataframe.
    For example, if position 1 repeats 2 times, position 2 repeats once, and position 3 repeats 2 times,
    the dictionary will be: {1: 1, 2: 2}.
    Positions are the numbers in the "position" column.
    :param df:
    :return: dictionary of number of repetitions
    """
    count_dict = {}
    for index, row in df.iterrows():
        # get the position
        position = row['position']
        # count how many times this position appears in the dataframe
        count = df['position'].value_counts()[position]
        # add the count to the dictionary
        if count in count_dict:
            count_dict[count] += 1
        else:
            count_dict[count] = 1
    return count_dict


def create_position_repetitions_graph(data: dict, gene: str, path=None) -> None:
    """
    Usage:
    gene = "SLC26A4"
    path_for_output = f"{PATH_TO_INBAR}results/position_repetition_for_gene_specific/"
    df = pd.read_csv(f"{PATH_TO_INBAR}results/gene_specific_df/{gene}_with_position.csv", header=0)
    position_repetition_dict = count_position_repetitions(df)
    create_position_repetitions_graph(position_repetition_dict, gene, path_for_output)
    """
    # Extract the keys (index) and values from your dictionary
    index = list(data.keys())
    values = list(data.values())
    # Create a bar graph
    plt.bar(index, values)
    # Add labels to the axes
    plt.xlabel('Index')
    plt.ylabel('Values')
    # Annotate the bars with their respective values
    for i, value in enumerate(values):
        plt.text(index[i], value, str(value), ha='center', va='bottom')
    # Add a title
    plt.title(f"Position repetitions for {gene}")
    # Show the bar graph
    # plt.show()
    plt.savefig(f"{path}position_repetition_graph_{gene}.png")


def create_list_of_variants_from_predictions_vs_real_csv(path):
    """This function creates a list of variants from the predictions_vs_real csv file."""
    df = pd.read_csv(f"{path}/predictions_vs_real_with_variant.csv", header=0)
    variant_list = df['variant'].tolist()
    return variant_list


def create_fasta_file_from_list_of_variants(variant_list, path):
    """This function creates a fasta file from a list of variants.
    The fasta file is saved in the path folder.
    FASTA format:
    >gene_HUMAN variant1 variant2 ...
    Args:
        variant_list: list of variants.
        path: the path to the folder where the fasta file will be saved.
    """
    gene = path.split("/")[-1]
    # create the fasta file/s, each file contains maximum 100 variants
    variants_counter = 0
    # get the sequence of the gene from uniport
    seq = uni.get_sequence(gene)
    while variants_counter < len(variant_list):
        variants = " ".join(variant_list[variants_counter:variants_counter+100])
        with open(f"{path}/{gene}_{variants_counter}.fasta", "w") as f:
            f.write(f">{gene}_HUMAN {variants}")
            f.write(f"\n{seq}")
        variants_counter += 100


def add_source_column_according_to_input(df: pd.DataFrame, source: str):
    """
    This function adds a column to the dataframe with the source of the data.
    :param df: dataframe with the data
    :param source: the source of the data
    :return: dataframe with the source column
    """
    df['source'] = source
    return df


if __name__ == "__main__":
    # path_to_folders = "/home/inbar/variants/clinvar_variants/"
    # df = pd.read_csv(f"/home/inbar/clinvar_variants_for_validation_set_without_redundancy1.csv", header=0)
    # write_to_folders.create_variant_folders_from_dataframe(path_to_folders, df)

    # # Get AF structure
    # copy_pdb_files_to_all_variant_folders_in_all_gene_folders("/home/inbar/variants/clinvar_variants/")

    # Run FoldX
    path = "/home/inbar/variants/clinvar_variants/"
    create_mutation_file_to_all_variants(path)

    # # Create fasta file for each gene
    # path = "predictions_vs_real"
    # genes = ["WFS1", "SLC26A4", "FGFR1", "MYO7A", "COL2A1", "COL4A5", "COL4A3"]
    # for gene in genes:
    #     path_for_gene = f"{path}/{gene}"
    #     variant_list = create_list_of_variants_from_predictions_vs_real_csv(path_for_gene)
    #     create_fasta_file_from_list_of_variants(variant_list, path_for_gene)

    # features_df = pd.read_csv(f"{PATH_TO_OUTPUT_FOLDER}features_for_GJB2.csv", header=0)
    # feat_ext_col_by_col.main(features_df, PATH_TO_OUTPUT_FOLDER, "GJB2_with_positions")

    # features_df = pd.read_csv(f"{PATH_TO_OUTPUT_FOLDER}features.csv", header=0)
    # features_df = create_sub_dataframe_according_to_gene(features_df, "GJB2")
    # # save as csv
    # features_df.to_csv(f"{PATH_TO_OUTPUT_FOLDER}GJB2_dvd_with_source.csv", index=False)
    #
    # dvd_df = pd.read_csv(f"{PATH_TO_OUTPUT_FOLDER}GJB2_dvd_with_source.csv", header=0)
    # feat_ext_col_by_col.main(dvd_df, PATH_TO_OUTPUT_FOLDER, "GJB2_dvd_with_patmut_features.csv")

    # #combine the dataframes of dvd and patmut
    # dvd_df = pd.read_csv(f"{PATH_TO_OUTPUT_FOLDER}GJB2_dvd_with_patmut_features.csv", header=0)
    # patmut_df = pd.read_csv(f"{PATH_TO_OUTPUT_FOLDER}GJB2_patmut_with_source.csv", header=0)
    # combined_df = pd.concat([dvd_df, patmut_df])
    # combined_df.to_csv(f"{PATH_TO_OUTPUT_FOLDER}GJB2_combined_with_source.csv", index=False)
    # #
    # add position column
    # combined_df = pd.read_csv(f"{PATH_TO_OUTPUT_FOLDER}gene_specific_df/GJB2_combined_with_source.csv", header=0)
    # combined_df = add_position_column(combined_df)
    # combined_df.to_csv(f"{PATH_TO_OUTPUT_FOLDER}gene_specific_df/GJB2_combined_with_source_and_position.csv", index=False)


    # # save the dictionary to a csv file in the path folder (from path variable)
    # position_repetition_df = pd.DataFrame.from_dict(position_repetition_dict, orient='index', columns=['count'])
    # position_repetition_df.to_csv(f"{path_for_output}/position_repetition_{gene}.csv", index=True)

    # features_df = pd.read_csv(f"{PATH_TO_NOVEL_BENIGN_VARIANTS}/Benign_for_gene_specific/combined.csv", header=0)
    # feat_ext_col_by_col.main(features_df, PATH_TO_NOVEL_BENIGN_VARIANTS)
    #
    # features_df = pd.read_csv(f"{PATH_TO_OUTPUT_FOLDER}features.csv", header=0)
    # feat_ext_col_by_col.main(features_df, PATH_TO_OUTPUT_FOLDER)
    # combined_df = pd.read_csv(f"/home/inbar/variants/Benign_for_gene_specific/combined.csv", header=0)
    # features_df = pd.read_csv(f"{PATH_TO_OUTPUT_FOLDER}features.csv", header=0)
    # benign_df = pd.read_csv(f"/home/inbar/variants/Benign_for_gene_specific/benign.csv", header=0)
    #
    # df = add_source_column(combined_df, features_df, benign_df)
    #
    # # write the dataframe to a csv file
    # df.to_csv(f"{PATH_TO_OUTPUT_FOLDER}combined_with_source.csv", index=False)

    # # save the dataframe to a csv file
    # combined_df.to_csv(f"{PATH_TO_NOVEL_BENIGN_VARIANTS}/Benign_for_gene_specific/combined.csv", index=False)

    # create_mutation_file_to_all_variants(f"{PATH_TO_VARIANTS_FOLDER}Benign_for_gene_specific/")

    # run_on_gene_folders(f"{PATH_TO_NOVEL_BENIGN_VARIANTS}/Benign_for_gene_specific/", ["MYO7A"])

    # for gene in ["COL2A1", "COL4A3", "COL4A5", "WFS1", "FGFR1", "SLC26A4", "MYO7A"]:
    #     df_with_position = add_position_column(pd.read_csv(f"{PATH_TO_OUTPUT_FOLDER}/gene_specific_df/{gene}.csv", header=0))
    #     df_with_position.to_csv(f"{PATH_TO_OUTPUT_FOLDER}/gene_specific_df/{gene}_with_position.csv", index=False)
    #
    #
    # copy_pdb_files_to_all_variant_folders(f"{PATH_TO_VARIANTS_FOLDER}Benign_for_gene_specific/GJB2/", "P29033", "GJB2")
    # create_foldx_mutation_to_variants_in_gene_folder(f"{PATH_TO_VARIANTS_FOLDER}Benign_for_gene_specific/GJB2/")
    # Extract all paths in the folder

    # paths = [d for d in os.listdir(f"{PATH_TO_VARIANTS_FOLDER}Benign_for_gene_specific/GJB2/") if os.path.isdir(f"{PATH_TO_VARIANTS_FOLDER}Benign_for_gene_specific/GJB2/{d}")]
    # for path in paths:
    #     run_command_line_programs.run_opra_oda_sasa(f"{PATH_TO_VARIANTS_FOLDER}Benign_for_gene_specific/GJB2/{path}")