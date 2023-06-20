import uniprot_info as uni
import add_mutation as add_mut
import dvd_data as dvd
import os
import subprocess as sp
import extract_features as ext_feat
import logging as log
import pandas as pd

PATH_TO_DATA_FOLDER = "/home/inbar/DVDdata/"
PATH_TO_VARIANTS_FOLDER = "/home/inbar/variants/"
PATHOGENICITY_TYPES = ["Benign", "Pathogenic"]


def main():
    """Main function of automation script."""
    # Create log file
    log.basicConfig(filename='automation.log', level=log.DEBUG, format='%(asctime)s %(levelname)s %(message)s')
    errors = []

    for pathogenicity in PATHOGENICITY_TYPES:

        path = f"{PATH_TO_DATA_FOLDER}{pathogenicity}/"
        os.makedirs(f"{PATH_TO_VARIANTS_FOLDER}{pathogenicity}/")

        for root, dirs, files in os.walk(path):
            current_path = os.getcwd()

            for i, file_name in enumerate(files):
                # if checkpoint['gene'] is not None and checkpoint['variant'] is not None:
                #     if i < checkpoint['variant']:
                #         continue

                os.chdir(current_path)

                if file_name.endswith(".csv"):  # If the file is a csv file
                    file_path = os.path.join(root, file_name)  # Get the path to the file
                    # Calculation using other functions
                    gene_name = dvd.extract_gene_name(file_name)  # Get gene name
                    "For each file (gene), create folder in variants folder"
                    os.chdir(f"{os.getcwd()}/variants/")
                    gene_folder_path = f"{PATH_TO_VARIANTS_FOLDER}{pathogenicity}/{gene_name}"
                    os.makedirs(gene_folder_path)  # Create folder for gene
                    os.chdir(f"{gene_folder_path}/")

                    "For each file (gene), create list of all mutations for that gene"
                    amino_acid_changes_list = dvd.make_list_of_variants_per_gene(
                        f"{PATH_TO_DATA_FOLDER}{pathogenicity}/dvd.v9.vtable.{gene_name}.csv")  # Create list of vars for that gene
                    amino_acid_changes_1_letter_list = dvd.convert_variants_to_1_letter(
                        amino_acid_changes_list)  # Convert to 1 letter code
                    print(f"{gene_name} {amino_acid_changes_1_letter_list}")

                    "For each mutation in amino_acid_changes_list, create folder, and calculate features"
                    for variant in amino_acid_changes_1_letter_list:
                        print(f"Starting {gene_name} {variant}")
                        # Check if variant is in the uniprot sequence
                        if uni.get_sequence_length(gene_name) < dvd.extract_number(variant):
                            error_message = f"""Variant {variant} is not in the uniprot of {gene_name} (seq too short)."""
                            errors.append(error_message)
                            print(error_message)
                            continue
                        if add_mut.check_aa_in_sequence(variant[0], uni.get_sequence(gene_name), dvd.extract_number(variant)):
                            pass
                        else:
                            error_message = f"Variant {variant} is not in the uniprot of {gene_name}."
                            errors.append(error_message)
                            print(error_message)
                            continue
                        # Get uniprot data and make sure aa in chain
                        uni_id = uni.get_uniprot_id(gene_name)
                        # Open folder for variant
                        variant_folder_name = f"{gene_name}_{uni_id}_{variant}"
                        var_path = add_mut.create_path_to_file(os.getcwd(),
                                        variant_folder_name)  # creates a path to the variant folder, and creates it in the computer
                        if var_path == "Error":
                            print("Error in creating path to variant folder, folder already exists. Skipping variant.")
                            continue
                        os.chdir(f"{var_path}/")
                        # Download PDB file from alpha fold
                        pdb_file = add_mut.get_structure_af(uni_id)

                        # Creates foldx mutant
                        add_mut.create_mutation_using_foldx(variant, uni_id)

                        # """Extract features"""
                        # # Extract plddt value
                        # plddt_residue_value = ext_feat.get_plddt(residue_num, pdb_file_path)  # Get plddt value for residue
                        #
                        # # ODA
                        # oda_wt = ext_feat.get_oda(variant, uni_id, "wt")
                        # oda_mut = ext_feat.get_oda(variant, uni_id, "mut")
                        # delta_oda = ext_feat.get_delta_oda(variant, uni_id)
                        #
                        # # OPRA
                        # opra_wt = ext_feat.get_opra(variant, uni_id, "wt")
                        # opra_mut = ext_feat.get_opra(variant, uni_id, "mut")
                        # delta_opra = ext_feat.get_delta_opra(variant, uni_id)
                        #
                        # # Hydrogen bonds
                        # hbonds_wt = ext_feat.get_hbonds(variant, uni_id, "wt")
                        # hbonds_mut = ext_feat.get_hbonds(variant, uni_id, "mut")
                        # delta_hbonds = ext_feat.get_delta_hbonds(variant, uni_id)
                        #
                        # # Substitution matrix value
                        # sub_value = ext_feat.get_sub_mat(aa1,aa2)
                        #
                        # # Save features as csv row
                        # features = {"Gene_name": gene_name, "Variant": variant, "UniProt_ID": uni_id,
                        #             "PLDDT": plddt_residue_value, "ODA_wt": oda_wt, "ODA_mut": oda_mut,
                        #             "Delta_ODA": delta_oda, "OPRA_wt": opra_wt, "OPRA_mut": opra_mut,
                        #             "Delta_OPRA": delta_opra, "HBonds_wt": hbonds_wt, "HBonds_mut": hbonds_mut,
                        #             "Delta_HBonds": delta_hbonds, "Substitution_matrix": sub_value}
                        # df_variant_row = save_data_as_df(pathogenicity, features)
                        # df = df.append(df_variant_row)

                        # Go back to gene folder
                        os.chdir(f"{gene_folder_path}/")  # This should be the last line in the loop

                        print(f"Finished {gene_name} {variant}")


    # print(f"""---------------------------------
    #     All the variants that are not in the uniprot sequence, and needs to be dealt with manually, are:
    #     {errors}
    #     ---------------------------------""")


def extract_features():
    pass


# def save_data_as_df(pathogenicity: str, features: dict):
#     """Save data of variant as df single row, with pathogenicity (Benign/pathogenic), and features.
#     args:
#         pathogenicity (str): Benign or Pathogenic
#         features (dict): Dictionary of features
#     returns:
#         df (pd.DataFrame): Dataframe of variant data (single row)
#         """
#     df = pd.DataFrame(features, index=[0])
#     df['Pathogenicity'] = pathogenicity


if __name__ == '__main__':
    main()