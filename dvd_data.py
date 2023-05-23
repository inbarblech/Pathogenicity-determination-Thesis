import csv
import os
import re

AMINO_ACIDS = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Glu': 'E', 'Gln': 'Q', 'Gly': 'G',
               'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S',
               'Trp': 'W', 'Tyr': 'Y', 'Val': 'V', 'Thr': 'T'}


def make_list_of_variants_per_gene(file_name):
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


def convert_variants_to_1_letter(variants):
    """Converts the given list of variants to 1 letter amino acid code."""
    converted_variants = []
    for var in variants:
        step_one = var[3:]
        step_two = step_one[:-3]
        one_letters_variant = f"{AMINO_ACIDS[var[:3]]}{step_two}{AMINO_ACIDS[step_one[-3:]]}"
        converted_variants.append(one_letters_variant)
    return converted_variants


def extract_number(string):
    return int(re.search(r'\d+', string).group())
