import pandas as pd
from data_retrievel_and_feature_extraction import uniprot_info as uni

# Function that creates a fasta file.
# FASTA format:
# >gene_HUMAN variant1 variant2 ...


def create_list_of_variants_from_predictions_vs_real_csv(path):
    """This function creates a list of variants from the predictions_vs_real csv file."""
    df = pd.read_csv(f"{path}")
    variant_list = df['variant'].tolist()
    return variant_list


def create_fasta_file_from_list_of_variants(variant_list, path, seq=None, gene=None):
    """This function creates a fasta file from a list of variants.
    The fasta file is saved in the path folder.
    FASTA format:
    >gene_HUMAN variant1 variant2 ...
    Args:
        variant_list: list of variants.
        path: the path to the folder where the fasta file will be saved.
        seq: the sequence of the gene. If None, the sequence will be retrieved from UniProt.
    """
    if seq is None or gene is None:
        gene = path.split("\\")[-1]
        seq = uni.get_sequence(gene)
    # create the fasta file/s, each file contains maximum 100 variants
    variants_counter = 0
    while variants_counter < len(variant_list):
        variants = " ".join(variant_list[variants_counter:variants_counter+100])
        with open(f"{path}\\{gene}_{variants_counter}.fasta", "w") as f:
            f.write(f">{gene}_HUMAN {variants}")
            f.write(f"\n{seq}")
        variants_counter += 100


if __name__ == "__main__":
    genes = ["GJB2"]
    file_with_variants_name = None  # Initiated inside the for loop
    output_path = f"C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\benchmarking\\MutPred2\\fasta_files"
    for gene in genes:
        seq = uni.get_sequence(gene)
        file_with_variants_name = f"{gene}_predictions_LOPO_XGB.csv"
        variants_list = create_list_of_variants_from_predictions_vs_real_csv(f"C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\predictions_vs_real\\{gene}\\{file_with_variants_name}")
        create_fasta_file_from_list_of_variants(variants_list, f"{output_path}\\", seq, gene)
