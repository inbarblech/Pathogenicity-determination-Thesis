import os

import pandas as pd
from biopandas.pdb import PandasPdb
from data_retrievel_and_feature_extraction import uniprot_info as uni

PATH_TO_VARIANTS_FOLDER = "/home/inbar/variants/"

AMINO_ACIDS = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Glu': 'E', 'Gln': 'Q', 'Gly': 'G',
               'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S',
               'Trp': 'W', 'Tyr': 'Y', 'Val': 'V', 'Thr': 'T'}


def convert_hgvs_protein_to_variant(hgvs_protein: str) -> str:
    """Convert the given hgvs protein to a variant.
    Example: hgvs_protein = 'p.Ala204Gly' -> variant = 'A204G'
    use the AMINO_ACIDS dictionary to convert the amino acids to 1 letter code."""
    aa1 = hgvs_protein[2:5]
    aa2 = hgvs_protein[-3:]
    variant = f"{AMINO_ACIDS[aa1]}{hgvs_protein[5:-3]}{AMINO_ACIDS[aa2]}"
    return variant


def convert_variants_list_to_1_letter(variants):
    """Converts the given list of variants to 1 letter amino acid code."""
    converted_variants = []
    for var in variants:
        one_letters_variant = convert_variant_to_1_letter(var)
        converted_variants.append(one_letters_variant)
    return converted_variants


def convert_1_letter_aa_to_3_letter(letter_variant):
    """Convert a single variant to 3 letter amino acid code.
    for example: A -> Ala"""
    # if letter_variant in AMINO_ACIDS.values(), return the key of the value
    three_letter_variant = [key for key, value in AMINO_ACIDS.items() if value == letter_variant][0]
    return three_letter_variant


def convert_variant_to_1_letter(variant):
    """Convert a single variant to 1 letter amino acid code."""
    step_one = variant[3:]
    step_two = step_one[:-3]
    one_letter_variant = f"{AMINO_ACIDS[variant[:3]]}{step_two}{AMINO_ACIDS[step_one[-3:]]}"
    return one_letter_variant


def get_aa1_aa2(variant: str) -> tuple:
    """Get the aa1 and aa2 from the variant
    Example: variant = 'A204G' -> aa1 = 'A', aa2 = 'G'"""
    aa1 = variant[0]
    aa2 = variant[-1]
    return aa1, aa2


def get_residue_number_from_variant_folder(folder_name: str) -> int:
    """Get the residue number from the folder name.
    Example: folder_name = 'ACTB_P60709_A204G' -> residue_number = 204"""
    residue_number = int(folder_name.split('_')[-1][1:-1])
    return residue_number


def get_variant_from_variant_folder(folder_name: str) -> str:
    """Get the variant from the folder name.
    Example: folder_name = 'ACTB_P60709_A204G' -> variant = 'A204G'"""
    variant = folder_name.split('_')[-1]
    return variant


def get_gene_from_variant_folder(folder_name: str) -> str:
    """Get the gene from the folder name.
    Example: folder_name = 'ACTB_P60709_A204G' -> gene = 'ACTB'"""
    gene = folder_name.split('_')[0]
    return gene


def get_amino_acid_of_variant_from_variant_folder(folder_name: str) -> str:
    """Get the amino acid of the variant from the folder name.
    Example: folder_name = 'ACTB_P60709_A204G' -> amino_acid = 'G'"""
    amino_acid = folder_name.split('_')[-1][-1]
    return amino_acid


def get_amino_acid_of_wt_from_variant_folder(folder_name: str) -> str:
    """Get the amino acid of the wt from the folder name.
    Example: folder_name = 'ACTB_P60709_A204G' -> amino_acid = 'A'"""
    amino_acid = folder_name.split('_')[-1][0]
    return amino_acid


def get_variant_location_from_variant_folder(folder_name: str) -> str:
    """Get the variant location from the folder name.
    Example: folder_name = 'ACTB_P60709_A204G' -> variant_location = '204'"""
    variant_location = folder_name.split('_')[-1][1:-1]
    return variant_location


def get_folder_name_from_path(path: str) -> str:
    """Get the folder name from the path.
    Example: path = '/home/inbar/variants/ACTB_P60709_A204G' -> folder_name = 'ACTB_P60709_A204G'"""
    path = path.strip()
    folder_name = os.path.basename(path.rstrip('/')).split('/')[-1]
    return folder_name


def get_folder_name_from_path_windows(path: str) -> str:
    path = path.strip()
    folder_name = os.path.basename(path.rstrip('\\')).split('\\')[-1]
    return folder_name


def extract_files_from_type(files_list: list, type: str) -> list:
    """Extract files from the given type from the given files list."""
    extracted_files = []
    for file in files_list:
        if file.endswith(type):
            extracted_files.append(file)
    return extracted_files


def convert_pdb_to_dataframe(pdb_path: str) -> pd.DataFrame:
    """Convert the given pdb file to a dataframe."""
    ppdb_df = PandasPdb().read_pdb(pdb_path).df['ATOM']
    return ppdb_df


def dataframe_to_pdb(dataframe, output_pdb_path):
    with open(output_pdb_path, 'w') as pdb_file:
        for index, row in dataframe.iterrows():
            pdb_line = f"ATOM  {index + 1:5} {row['NAME']:^4} {row['RESNAME']:^3} {row['CHAINID']}{row['RESIDUESEQ']:4}    {row['X']:8.3f}{row['Y']:8.3f}{row['Z']:8.3f}{row['OCCUPANCY']:6.2f}{row['TEMPFACTOR']:6.2f}          {row['ELEMENT']:>2}\n"
            pdb_file.write(pdb_line)


def get_variant_path(gene, variant, pathogenicity):
    """Get the variant path from the given gene, variant and pathogenicity."""
    # convert variant from E117K to ACTB_P60709_E117K by adding the gene name and the uniprot id
    uniprot_id = uni.get_uniprot_id(gene)
    variant_path = f"{gene}_{uniprot_id}_{variant}"
    return f"{PATH_TO_VARIANTS_FOLDER}{pathogenicity}/{gene}/{variant_path}/"


def get_uniprot_id_from_path(path: str) -> str:
    """Get the uniprot id from the given path."""
    variant_folder = path.split('/')[-2]
    uniprot_id = variant_folder.split('_')[1]
    return uniprot_id


def get_gene_name_from_path(path: str) -> str:
    """Get the gene name from the given path."""
    variant_folder = path.split('/')[-2]
    gene_name = variant_folder.split('_')[0]
    return gene_name


def get_uniprot_ids_from_df(df: pd.DataFrame) -> list:
    genes_set = get_genes_set_from_df(df)
    uniprot_ids = []
    for gene in genes_set:
        uni_id = uni.get_uniprot_id(gene)
        uniprot_ids.append(uni_id)
    return uniprot_ids


def get_uniprot_ids_and_gene_names_from_df(df: pd.DataFrame) -> pd.DataFrame:
    genes_set = get_genes_set_from_df(df)
    uniprot_ids = pd.DataFrame(columns=['gene', 'uniprot_id'])
    for gene in genes_set:
        uni_id = uni.get_uniprot_id(gene)
        uniprot_ids = uniprot_ids.append({'gene': gene, 'uniprot_id': uni_id}, ignore_index=True)
    return uniprot_ids


def get_genes_set_from_df(df: pd.DataFrame) -> set:
    genes = df['gene']
    genes_set = set(genes)
    return genes_set


def get_number_of_variants_per_gene_dict(df: pd.DataFrame) -> dict:
    """Recieves a dataframe with one type of pathogenicity.
    Returns dictionary of number of variants per gene"""
    gene_set = get_genes_set_from_df(df)
    number_of_variants_per_gene = dict()
    for gene in gene_set:
        number_of_variants = len(df[df["gene"] == gene])
        number_of_variants_per_gene[gene] = number_of_variants
    return number_of_variants_per_gene


def split_dataframe_by_group(dataframe, group_column, group_value1, group_value2):
    """
    Split a DataFrame into two sub-DataFrames based on the values in the 'group_column'.

    Parameters:
    - dataframe: DataFrame to be split.
    - group_column: The name of the column used for splitting.
    - group_value1: The value in 'group_column' to include in the first sub-DataFrame.
    - group_value2: The value in 'group_column' to include in the second sub-DataFrame.

    Returns:
    - Two sub-DataFrames.
    """
    df_group1 = dataframe[dataframe[group_column] == group_value1]
    df_group2 = dataframe[dataframe[group_column] == group_value2]
    return df_group1, df_group2


def get_gene_path_from_variant_path(variant_path: str) -> str:
    """Get the gene path from the given variant path, for example:
    variant_path = '/home/inbar/variants/pathogenic/ACTB/ACTB_P60709_A204G'
    gene_path = '/home/inbar/variants/pathogenic/ACTB'"""
    gene_path = variant_path[:variant_path.rfind('/')]
    return gene_path
