

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

