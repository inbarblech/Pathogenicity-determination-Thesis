import os
import uniprot_info as uni
import pandas as pd
import add_mutation as add_mut


def create_variant_folders_from_dataframe(path_to_genes_folder, variants_df):
    """Create variant folders in the path from the data in the given variants dataframe.
    The dataframe should contain the following columns:
    - gene
    - variant
    - pathogenicity
    - uniprot_id
    For example, for the variant E117K in ACTB, the path should be:
    /path_to_genes_folder/ACTB/ACTB_P60709_E117K/
    Assumes that the gene folder ('ACTB') already exists.
    """
    for index, row in variants_df.iterrows():
        gene = row['gene']
        variant = row['variant']
        uniprot_id = row['uniprot_id']
        variant_path = f"{gene}_{uniprot_id}_{variant}"
        variant_folder_path = f"{path_to_genes_folder}{gene}/{variant_path}/"
        if not os.path.exists(variant_folder_path):
            os.makedirs(variant_folder_path)
        else:
            print(f"Variant folder already exists: {variant_folder_path}")


def copy_af_file_to_all_subfolders(path_to_gene_folder):
    """This function copies the .pdb file in the path to all the subfolders of the given path."""
    for root, dirs, files in os.walk(path_to_gene_folder):
        for file in files:
            if file.endswith(".pdb"):
                pdb_file_path = os.path.join(root, file)
                for root2, dirs2, files2 in os.walk(root):
                    for dir2 in dirs2:
                        if not dir2.endswith("PDBs"):
                            new_pdb_file_path = os.path.join(root2, dir2, file)
                            print(f"Copying {pdb_file_path} to {new_pdb_file_path}")
                            os.system(f"cp {pdb_file_path} {new_pdb_file_path}")


def check_if_aa_is_in_af(data: pd.DataFrame) -> None:
    """Checks if the amino acid in the variant is in the sequence of the protein from AlphaFold."""
    errors = pd.DataFrame(columns=['gene', 'variant'])
    # create a dictionary with the gene name as the key and the variants as the value.
    gene_variants_dict = {}
    #extract the length of data
    length = len(data)
    print(f"Going over {length} variants")
    count = 0
    for index, row in data.iterrows():
        gene = row['gene']
        variant = row['variant']
        if gene not in gene_variants_dict:
            gene_variants_dict[gene] = [variant]
        else:
            gene_variants_dict[gene].append(variant)
    genes = list(gene_variants_dict.keys())
    # for each gene, check if the variant is in the sequence
    for gene in genes:
        uniprot_id = uni.get_uniprot_id(gene)
        for variant in gene_variants_dict[gene]:
            variant_aa = variant[0]  # get the amino acid from the variant, e.g. "A" from "A123B"
            pos = int(variant[1:-1])  # get the position from the variant, e.g. 123 from "A123B"
            seq = uni.get_sequence(gene)
            count += 1
            print(f"Checking variant {count} out of {length}")
            if uni.get_sequence_length(gene) < int(pos):  # check if the position is in the sequence
                error = "Variant exceeds sequence length"
                # add variant to errors dataframe
                errors = errors.append({'gene': gene, 'variant': variant, 'uniprot_id': uniprot_id,
                                        'error': error}, ignore_index=True)
                continue
            if add_mut.check_aa_in_sequence(variant_aa, seq, pos) is False:
                error = "Variant not in sequence"
                errors = errors.append({'gene': gene, 'variant': variant, 'uniprot_id': uniprot_id,
                                        'error': error}, ignore_index=True)
                continue
    # write errors to file in gene folder
    errors.to_csv("/home/inbar/variants/Benign_for_gene_specific/not_in_sequence_error.csv")


if __name__ == "__main__":
    path = "/home/inbar/variants/Benign_for_gene_specific/"
    # benign_df = pd.read_csv("/home/inbar/variants/Benign_for_gene_specific/benign.csv")
    # create_variant_folders_from_dataframe(path, benign_df)

    # create a list of all gene folders paths inside the given path
    # gene_folders_paths = []
    # for root, dirs, files in os.walk(path):
    #     for dir in dirs:
    #         gene_folders_paths.append(os.path.join(root, dir))
    # for path in gene_folders_paths:
    #     copy_af_file_to_all_subfolders(path)

    benign_df = pd.read_csv("/home/inbar/variants/Benign_for_gene_specific/benign.csv")
    check_if_aa_is_in_af(benign_df)
