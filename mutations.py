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

############################################## Search UniProt ID's without MSA file ###################################
dir_to_msas = "C:\\Users\\InbarBlech\\Downloads\\arffs"
uniprot_id_list = ['Q6ISB3', 'Q8N205', 'Q15031', 'Q7RTW8', 'Q14055', 'O75616', 'P58418', 'P21281', 'P35240', 'Q9HAB3',
                   'Q7Z412', 'O95677', 'O00462', 'P12107', 'P56693', 'P02458', 'Q14123', 'Q9NPC8', 'P22607', 'Q14651',
                   'Q96QU1', 'Q9Y4X0', 'Q9ULP9', 'Q14031', 'Q01668', 'P11487', 'Q8TF64', 'Q53GD3', 'Q9NQ40', 'Q9H2D6',
                   'O15160', 'P13942', 'O95484', 'O00400', 'Q8WZ04', 'Q8N2K0', 'Q9HCJ1', 'P0DPB6', 'O60220', 'Q9NSB8',
                   'P11362', 'Q86X52', 'Q969T9', 'Q12951', 'Q9NR28', 'Q8WZ55', 'Q00604', 'Q9UM54', 'O95136', 'Q01955',
                   'Q9HBG4', 'P29400', 'P51787', 'Q8TAF8', 'O95831', 'Q13428', 'Q2WEN9', 'O43435', 'Q08828', 'Q9UHP9',
                   'Q16740', 'Q9NZA1', 'O15041', 'P08581', 'Q9NSV4', 'O60779', 'Q8N196', 'Q7RTU9', 'Q4KMQ1', 'Q9Y6N9',
                   'P29033', 'Q9H9Y6', 'Q13608', 'O75712', 'Q8IU57', 'Q9UBR4', 'Q01814', 'Q9H015', 'P14210', 'Q8IVV2',
                   'P15313', 'P21583', 'Q15475', 'Q8N5K1', 'Q15319', 'Q7Z5J4', 'O75838', 'P63261', 'O43364', 'Q13402',
                   'P14653', 'P58743', 'Q9UDY2', 'Q96P20', 'P07196', 'Q14050', 'O43933', 'Q01973', 'Q8IVW8', 'P49335',
                   'P23771', 'Q13588', 'Q96FG2', 'O95500', 'Q14894', 'O76024', 'Q8N4S9', 'Q9NZW4', 'P22570', 'O00754',
                   'P56178', 'P35579', 'Q9H061', 'Q9C091', 'Q9Y276', 'P56696', 'P78508', 'A8MXD5', 'P15863', 'Q96RR1',
                   'Q9Y4F9', 'O60610', 'Q5JTW2', 'P35241', 'O43623', 'Q9H1P3', 'B1AK53', 'O43511', 'P51659', 'Q8NB90',
                   'P15382', 'P08493', 'O43405', 'P58215', 'P49590', 'Q8TCS8', 'Q9UBL9', 'P35237', 'Q15046', 'Q9UMZ3',
                   'P60709', 'Q96I59', 'Q9P202', 'Q86SU0', 'Q9P2R7', 'O75443', 'Q13253', 'Q9NRS4', 'P53420', 'O75030',
                   'P26358', 'Q12929', 'Q8TE12', 'P24821', 'Q8IVM0', 'Q8NBS3', 'Q8NEV4', 'Q8NDX2', 'P48740', 'Q8TDI8',
                   'Q5H9S7', 'Q13480', 'Q99502', 'P14138', 'Q04900', 'P60891', 'Q9H5Y7', 'O43314', 'Q9UHG0', 'Q7Z406',
                   'Q13127', 'O95718', 'O60443', 'P81274', 'P43251', 'Q9H6S3', 'P23760', 'Q495M9', 'Q99952', 'Q8NEW7',
                   'P20849', 'P10589', 'Q96D09', 'Q0ZLH3', 'Q9HC10', 'Q9UNH5', 'Q6IEE7', 'Q9H5P4', 'A6NFK2', 'Q8N6M3',
                   'O60487', 'P05549', 'O60907', 'P21802', 'O60313']
uniprot_list_without_msa = get_uniprots_without_msa(uniprot_id_list, dir_to_msas)

with open("C:\\Users\\InbarBlech\\Downloads\\uniprots_without_msa.txt", 'w') as f:
    for gene in uniprot_list_without_msa:
        f.write(f"{gene}\n")
def get_uniprots_without_msa(uniprot_id_list, dir):
    uniprots_without_msa = []
    uniprot_in_directory = find_msa_files(dir)
    for uniprot_id in uniprot_id_list:
        if uniprot_id not in uniprot_in_directory:
            uniprots_without_msa.append(uniprot_id)
    return uniprots_without_msa


def find_msa_files(directory):
    msa_files = []  # Create an empty list to store the file names
    # Walk through the directory and its subdirectories
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(".msa"):
                # Remove the ".msa" extension and append to the list
                file_name_without_extension = os.path.splitext(file)[0]
                msa_files.append(file_name_without_extension)
    return msa_files

######################################################################################################################