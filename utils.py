from data_retrievel_and_feature_extraction import uniprot_info as uni


def check_variant_valid(variant, gene):
    """
    Check if the variant is valid for the given gene.
    Args:
        variant (str): The variant to check.
        gene (str): The gene to check.
    Returns:
        bool: True if the variant is valid, False otherwise
    """
    # Check if the variant is in the correct format
    if not variant.isalnum():
        return False
    # Check if the amino acids are valid
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    if variant[0] not in amino_acids or variant[-1] not in amino_acids:
        return False
    # Check if the position is valid (all characters except the first and last)
    if not variant[1:-1].isdigit():
        return False

    # Check if this position is valid for the given gene, using uniprot
    gene_sequence = uni.get_sequence(gene)
    # Check if the position is within the length of the gene sequence
    if int(variant[1:-1]) > len(gene_sequence):
        return False
    # Check if the amino acid at this position is the same as the one in the variant
    if gene_sequence[int(variant[1:-1]) - 1] != variant[0]:
        return False

    return True



