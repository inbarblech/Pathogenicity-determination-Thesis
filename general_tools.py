
import os

AMINO_ACIDS = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Glu': 'E', 'Gln': 'Q', 'Gly': 'G',
               'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S',
               'Trp': 'W', 'Tyr': 'Y', 'Val': 'V', 'Thr': 'T'}


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
    three_letter_variant = AMINO_ACIDS[letter_variant]
    return three_letter_variant


def convert_variant_to_1_letter(variant):
    """Convert a single variant to 1 letter amino acid code."""
    step_one = variant[3:]
    step_two = step_one[:-3]
    one_letter_variant = f"{AMINO_ACIDS[variant[:3]]}{step_two}{AMINO_ACIDS[step_one[-3:]]}"
    return one_letter_variant


def get_aa1_aa2(variant: str) -> tuple:
    """Get the aa1 and aa2 from the variant"""
    aa1 = variant[0]
    aa2 = variant[-1]
    return aa1, aa2


def get_residue_number_from_variant_folder(folder_name: str) -> int:
    """Get the residue number from the folder name.
    Example: folder_name = 'ACTB_P60709_A204G' -> residue_number = 204"""
    residue_number = int(folder_name.split('_')[-1][1:-1])
    return residue_number


def get_amino_acid_of_variant_from_variant_folder(folder_name: str) -> str:
    """Get the amino acid of the variant from the folder name.
    Example: folder_name = 'ACTB_P60709_A204G' -> amino_acid = 'G'"""
    amino_acid = folder_name.split('_')[-1][-1]
    return amino_acid


def get_variant_location_from_variant_folder(folder_name: str) -> str:
    """Get the variant location from the folder name.
    Example: folder_name = 'ACTB_P60709_A204G' -> variant_location = '204'"""
    variant_location = folder_name.split('_')[-1][1:-1]
    return variant_location


def get_folder_name_from_path(path: str) -> str:
    """Get the folder name from the path."""
    folder_name = os.path.basename(path)
    return folder_name


def extract_files_from_type(files_list: list, type: str) -> list:
    """Extract files from the given type from the given files list."""
    extracted_files = []
    for file in files_list:
        if file.endswith(type):
            extracted_files.append(file)
    return extracted_files
