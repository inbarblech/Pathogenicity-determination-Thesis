# TRASH


# main.py
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



# def main():
#     """Main function of automation script."""
#     # Create log file
#     log.basicConfig(filename='automation.log', level=log.DEBUG, format='%(asctime)s %(levelname)s %(message)s')
#     errors = []
#
#     for pathogenicity in PATHOGENICITY_TYPES:
#         path = f"{PATH_TO_DATA_FOLDER}{pathogenicity}/"
#         os.makedirs(f"{PATH_TO_VARIANTS_FOLDER}{pathogenicity}/")
#
#         for root, dirs, files in os.walk(path):
#             current_path = os.getcwd()
#
#             for file_name in files:
#                 os.chdir(current_path)
#
#                 if file_name.endswith(".csv"):  # If the file is a csv file
#                     file_path = os.path.join(root, file_name)  # Get the path to the file
#                     # Calculation using other functions
#                     gene_name = dvd.extract_gene_name(file_name)  # Get gene name
#                     "For each file (gene), create folder in variants folder"
#                     os.chdir(f"{os.getcwd()}/variants/")
#                     gene_folder_path = f"{PATH_TO_VARIANTS_FOLDER}{pathogenicity}/{gene_name}"
#                     os.makedirs(gene_folder_path)  # Create folder for gene
#                     os.chdir(f"{gene_folder_path}/")
#
#                     "For each file (gene), create list of all mutations for that gene"
#                     amino_acid_changes_list = dvd.make_list_of_variants_per_gene(
#                         f"{PATH_TO_DATA_FOLDER}{pathogenicity}/dvd.v9.vtable.{gene_name}.csv")  # Create list of vars for that gene
#                     amino_acid_changes_1_letter_list = dvd.convert_variants_to_1_letter(
#                         amino_acid_changes_list)  # Convert to 1 letter code
#                     print(f"{gene_name} {amino_acid_changes_1_letter_list}")
#
#                     "For each mutation in amino_acid_changes_list, create folder, and calculate features"
                    for variant in amino_acid_changes_1_letter_list:
                        print(f"Starting {gene_name} {variant}")
                        # Check if variant is in the uniprot sequence
                        if uni.get_sequence_length(gene_name) < dvd.extract_number(variant):
                            error_message = f"""Variant {variant} is not in the uniprot of {gene_name} (seq too short)."""
                            errors.append(error_message)
                            print(error_message)
                            continue
                        if add_mut.check_aa_in_sequence(variant[0], uni.get_sequence(gene_name),
                                                        dvd.extract_number(variant)):
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

                        # Go back to gene folder
                        os.chdir(f"{gene_folder_path}/")  # This should be the last line in the loop

                        print(f"Finished {gene_name} {variant}")

    # Print errors









#
# def get_plddt(residue_num: str, af_file_path: str):
#     """returns the pLDDT score of the residue
#     Uses the alphafold pdb file to extract the pLDDT score"""
#
#     # convert file to csv
#     af_df = tools.convert_pdb_to_dataframe(af_file_path)
#
#     # get the pLDDT score, which is the first value in the 11th column where the residue number = residue_num
#     # in the dataframe
#     value = af_df[af_df[5] == residue_num][10].iloc[0]
#
#     # value = None  # The value to be returned
#     # # read the csv file
#     # with open(pdb_file_pwd, 'r') as csv_file:
#     #     reader = csv.reader(csv_file)
#     #     for row in reader:
#     #         if row[5] == residue_num:
#     #             value = row[10]
#     #             break  # Exit the loop after finding the first row
#     # if value is None:
#     #     raise ValueError("The residue number is not valid, or the residue is not in the pdb file.")
#
#     return value


