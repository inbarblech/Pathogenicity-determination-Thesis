import uniprot_info as uni
import add_mutation as add_mut
import dvd_data as dvd
import os
import subprocess as sp
import extract_features as ext_feat


def main():
    """Main function of automation script."""
    # Go through the benign genes folder
    errors = []
    for pathogenicity in ["Benign", "Pathogenic"]:
        path = f"/home/inbar/DVDdata/{pathogenicity}/"
        os.makedirs(f"/home/inbar/variants/{pathogenicity}/")
        for root, dirs, files in os.walk(path):
            current_path = os.getcwd()
            for file_name in files:
                os.chdir(current_path)
                if file_name.endswith(".csv"):
                    file_path = os.path.join(root, file_name)
                    # Calculation using other functions
                    gene_name = dvd.extract_gene_name(file_name)  # Get gene name
                    "For each file (gene), create folder in variants folder"
                    os.chdir(f"{os.getcwd()}/variants/")
                    gene_folder_path = f"/home/inbar/variants/{pathogenicity}/{gene_name}"
                    os.makedirs(gene_folder_path)  # Create folder for gene
                    os.chdir(f"{gene_folder_path}/")

                    "For each file (gene), create list of all mutations for that gene"
                    amino_acid_changes_list = dvd.make_list_of_variants_per_gene(
                        f"/home/inbar/DVDdata/{pathogenicity}/dvd.v9.vtable.{gene_name}.csv")  # Create list of vars for that gene
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
    
                        # Create feature files and extract features
                        add_mut.create_mutation_using_foldx(variant, uni_id)
                        ext_feat.get_oda(uni_id)
                        os.chdir(f"{gene_folder_path}/")  # This should be the last line in the loop
                        # Make calculations (Features) including deltas, and save as series/dictionary.

                        # Append row to csv/df, with info as pathogenicity (Benign/pathogenic), and features.
    print(f"""
        ---------------------------------
        All the variants that are not in the uniprot sequence, and needs to be dealt with manually, are:
        {errors}
        ---------------------------------""")

if __name__ == '__main__':
    main()
