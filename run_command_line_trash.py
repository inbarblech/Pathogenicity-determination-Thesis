def run_oda(variant_path: str):
    """Creates a file with the ODA score for the given PDB file.
    The file is created in the variant folder.
    Args:
        variant_path (str): The path to the variant folder.
    """
    # change the working directory to the variant folder
    os.chdir(variant_path)

    # get gene_id from folder name
    folder_name = tools.get_folder_name_from_path(variant_path)
    variant_location = tools.get_variant_location_from_variant_folder(folder_name)
    gene_id = folder_name.split('_')[1]

    # Run oda for mut
    amino_acid = tools.get_amino_acid_of_variant_from_variant_folder(folder_name)
    three_letter_amino_acid = tools.convert_1_letter_aa_to_3_letter(amino_acid)
    # make three_letter_amino_acid uppercase
    three_letter_amino_acid = three_letter_amino_acid.upper()

    # sasa_script = "/home/inbar/sasa.sh"
    # run_pyDock_oda_mut_script = "/home/inbar/PyDock/PyDock3/run_pyDock_oda_mut.sh"
    # run_pyDock_oda_wt_script = "/home/inbar/PyDock/PyDock3/run_pyDock_oda_wt.sh"

    # run oda
    # command = f'bash -c "source {sasa_script} && run_pyDock_oda_mut \'{three_letter_amino_acid}\' \'{variant_location}\' \'{gene_id}\'"'
    #sp.run(command)
    # command = f'bash -c "source {sasa_script} && run_pyDock_oda_mut \'{gene_id}\'"'
    #sp.run(command)


def run_opra(variant_path: str):
    """creates a file with the OPRA score for the given PDB file."""
    # change the working directory to the variant folder
    os.chdir(variant_path)

    # get gene_id from folder name
    folder_name = tools.get_folder_name_from_path(variant_path)
    variant_location = tools.get_variant_location_from_variant_folder(folder_name)
    gene_id = folder_name.split('_')[1]

    # Run opra for wt and mutq
    amino_acid = tools.get_amino_acid_of_variant_from_variant_folder(folder_name)
    three_letter_amino_acid = tools.convert_1_letter_aa_to_3_letter(amino_acid)
    # make three_letter_amino_acid uppercase
    three_letter_amino_acid = three_letter_amino_acid.upper()
    # run opra
    command = f"run_pyDock_opra_mut {three_letter_amino_acid} {variant_location} {gene_id}"
    sp.run(command, shell=True)
    command = f"run_pyDock_opra_wt {gene_id}"
    sp.run(command, shell=True)

def run_sasa(path_to_variant_folder: str):
    """Creates a file with the SASA score for the given PDB file."""
    os.chdir(path_to_variant_folder)
    # get gene_id from folder name
    folder_name = tools.get_folder_name_from_path(path_to_variant_folder)
    variant_location = tools.get_variant_location_from_variant_folder(folder_name)
    gene_id = folder_name.split('_')[1]

    command = f"run_sasa_wt {gene_id}"
    sp.run(command, shell=True)
    amino_acid = tools.get_amino_acid_of_variant_from_variant_folder(folder_name)
    three_letter_amino_acid = tools.convert_1_letter_aa_to_3_letter(amino_acid)
    # make three_letter_amino_acid uppercase
    three_letter_amino_acid = three_letter_amino_acid.upper()
    # run sasa
    command = f"run_sasa_mut {three_letter_amino_acid} {variant_location} {gene_id}"
    sp.run(command, shell=True)