import csv
import os
import re

PATH_TO_DATA_FOLDER = "/home/inbar/DVDdata/"
PATH_TO_CSV_BENIGN = "/home/inbar/all_BENIGN_variants.csv"
PATH_TO_CSV_PATHOGENIC = "/home/inbar/all_PATHOGENIC_variants.csv"
PATH_TO_VARIANTS_FOLDER = "/home/inbar/variants/"
PATHOGENICITY_TYPES = ["Benign", "Pathogenic"]


def make_list_of_variants_for_all_genes_from_dvd_data(path_to_folder_pathogenicity):
    """Creates a list of variants from the given folder (all genes!).
    the variants are given under the column named 'HGVS Protein Change"""
    variants = []
    for root, dirs, files in os.walk(path_to_folder_pathogenicity):
        for file_name in files:
            if file_name.endswith(".csv"):
                file_path = os.path.join(root, file_name)
                variants.extend(make_list_of_variants_per_gene_from_dvd_data(file_path))
    with open('variants.txt', 'w') as file:
        for variant in variants:
            file.write(variant + '\n')
    return variants


def count_all_variants_in_folder(path_to_folder):
    """Counts the overall variants in the given folder separately for benign and pathogenic folders.
    USE: benign, pathogenic = count_all_variants_in_folder(PATH_TO_DATA_FOLDER)"""
    benign_variants = set()
    pathogenic_variants = set()

    for root, dirs, files in os.walk(path_to_folder):
        for file_name in files:
            if file_name.endswith(".csv"):
                file_path = os.path.join(root, file_name)
                with open(file_path, 'r') as file:
                    file.readline()  # Skip the first line
                    reader = csv.DictReader(file)
                    for row in reader:
                        variant = row['HGVS Protein Change']
                        if "Benign" in root:
                            benign_variants.add(variant)
                        elif "Pathogenic" in root:
                            pathogenic_variants.add(variant)

    return len(benign_variants), len(pathogenic_variants), benign_variants, pathogenic_variants


def make_list_of_variants_per_gene_from_dvd_data(file_name):
    """Creates a list of variants from the given csv file (per one gene!).
    the variants are given under the column named 'HGVS Protein Change"""
    variants = []
    aa_changes = []
    with open(file_name, 'r') as file:
        file.readline()
        reader = csv.DictReader(file)
        for row in reader:
            variants.append(row['HGVS Protein Change'])
    for var in variants:
        aa_change = var.split(':p.')[1]
        aa_changes.append(aa_change)
    return aa_changes


def extract_gene_name(file_name, prefix="dvd.v9.vtable."):
    if file_name.startswith(prefix):
        remaining = file_name[len(prefix):]
        dot_index = remaining.find('.')
        if dot_index != -1:
            gene_name = remaining[:dot_index]
            return gene_name
    return None


def extract_number(string):
    return int(re.search(r'\d+', string).group())


if __name__ == "__main__":
    benign_count, pathogenic_count, benign, pathogenic = count_all_variants_in_folder(PATH_TO_DATA_FOLDER)
    # pathogenic2 = count_all_variants_from_dvd_data_in_one_pathogenicity_folder(f"{PATH_TO_DATA_FOLDER}PATHOGENICITY_TYPES[1]")

    print(f"Number of benign variants: {benign_count}")
    print(f"Number of pathogenic variants: {pathogenic_count}")
    print(f"Number of variants: {benign_count + pathogenic_count}")
    print(f"Number of variants in both benign and pathogenic: {len(benign.intersection(pathogenic))}")
    print(f"Number of variants in benign but not in pathogenic: {len(benign.difference(pathogenic))}")
    print(f"Number of variants in pathogenic but not in benign: {len(pathogenic.difference(benign))}")
