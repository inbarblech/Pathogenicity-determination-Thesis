import csv
import os
import re


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