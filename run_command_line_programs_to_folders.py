
import subprocess as sp
import os
import general_tools as tools


def run_opra_oda_sasa(variant_path: str):
    """creates a file with the OPRA, ODA and SASA score for the given PDB files of wt and mut, if the files exists."""
    # change the working directory to the variant folder
    os.chdir(variant_path)

    # get gene_id from folder name
    folder_name = tools.get_folder_name_from_path(variant_path)
    variant_location = tools.get_variant_location_from_variant_folder(folder_name)
    gene_id = folder_name.split('_')[1]

    # Check if there is an AF file in this folder, then run all the functions for wt
    # if not os.path.isfile(f"AF_{gene_id}.pdb"):
    #     raise Exception(f"Error, no AF file for variant {variant_path}")
    # else:
    #     # Run oda for wt
    #     command_wt = f"pyDock3 AF_{gene_id}.pdb oda"
    #     sp.run(command_wt, shell=True)
    #     # Run opra for wt
    #     # command_wt = f"pyDock3 AF_{gene_id}.pdb opra"
    #     # sp.run(command_wt, shell=True)
    #     # Run sasa for wt
    #     command_wt = f"dr_sasa -m 0 -i AF_{gene_id}.pdb"
    #     sp.run(command_wt, shell=True)

    # Run for mut
    amino_acid = tools.get_amino_acid_of_variant_from_variant_folder(folder_name)
    three_letter_amino_acid = tools.convert_1_letter_aa_to_3_letter(amino_acid)
    # make three_letter_amino_acid uppercase
    three_letter_amino_acid = three_letter_amino_acid.upper()
    # Check if there is a mut file in this folder, then run all the functions for mut
    if not os.path.isfile(f"{three_letter_amino_acid}{variant_location}_AF_{gene_id}.pdb"):
        raise Exception(f"Error, no mutant pdb in {variant_path}")
    else:
        # Run oda for mut
        command_mut = f"pyDock3 {three_letter_amino_acid}{variant_location}_AF_{gene_id}.pdb oda"
        sp.run(command_mut, shell=True)
        # Run opra for mut
        # command_mut = f"pyDock3 {three_letter_amino_acid}{variant_location}_AF_{gene_id}.pdb opra"
        # sp.run(command_mut, shell=True)
        # Run sasa for mut
        command_mut = f"dr_sasa -m 0 -i {three_letter_amino_acid}{variant_location}_AF_{gene_id}.pdb"
        sp.run(command_mut, shell=True)


def run_sasa(variant_path: str, file_name: str, mut_file_name = None) -> str:
    """Creates a file with the sasa score for the given PDB file.
    The file is created in the variant folder.
    This function can run sasa for wt and mut, if both files are given.
    Args:
        variant_path (str): The path to the variant folder.
        file_name (str): The name of the PDB file to run sasa on.
        mut_file_name (str): The name of the PDB file to run sasa on of the mutation, if it exists.
    """
    # change the working directory to the variant folder
    os.chdir(variant_path)

    # Check if there is an AF file in this folder, then run oda for wt
    if not os.path.isfile(file_name):
        return f"Error, no AF file for variant {variant_path}"
    else:
        command_wt = f"dr_sasa -m 0 -i {file_name}"
        sp.run(command_wt, shell=True)

    if mut_file_name is not None:
        if not os.path.isfile(mut_file_name):
            return f"Error, no mut file for variant {variant_path}"
        else:
            command_mut = f"dr_sasa -m 0 -i {mut_file_name}"
            sp.run(command_mut, shell=True)

    return f"Successfully created SASA file for variant {variant_path}"


def run_oda(variant_path: str, file_name: str, mut_file_name = None) -> str:
    """Creates a file with the oda score for the given PDB file.
    The file is created in the variant folder.
    This function can run oda for wt and mut, if both files are given.
    Args:
        variant_path (str): The path to the variant folder.
        file_name (str): The name of the PDB file to run sasa on.
        mut_file_name (str): The name of the PDB file to run sasa on of the mutation, if it exists.
    """
    # change the working directory to the variant folder
    os.chdir(variant_path)

    # Check if there is an AF file in this folder, then run oda for wt
    if not os.path.isfile(file_name):
        return f"Error, no AF file for variant {variant_path}"
    else:
        command_wt = f"pyDock3 {file_name} oda"
        sp.run(command_wt, shell=True)

    if mut_file_name is not None:
        if not os.path.isfile(mut_file_name):
            return f"Error, no mut file for variant {variant_path}"
        else:
            command_mut = f"pyDock3 {mut_file_name} oda"
            sp.run(command_mut, shell=True)

    return f"Successfully created oda file for variant {variant_path}"


def run_opra(variant_path: str, file_name: str, mut_file_name=None) -> str:
    """Creates a file with the oda score for the given PDB file.
    The file is created in the variant folder.
    This function can run oda for wt and mut, if both files are given.
    Args:
        variant_path (str): The path to the variant folder.
        file_name (str): The name of the PDB file to run sasa on.
        mut_file_name (str): The name of the PDB file to run sasa on of the mutation, if it exists.
    """
    # change the working directory to the variant folder
    os.chdir(variant_path)

    # Check if there is an AF file in this folder, then run oda for wt
    if not os.path.isfile(file_name):
        return f"Error, no AF file for variant {variant_path}"
    elif os.path.isfile(f"{file_name}.opra"):
        return f"Function already run for variant {variant_path}"
    else:
        command_wt = f"pyDock3 {file_name} opra"
        sp.run(command_wt, shell=True)

    if mut_file_name is not None:
        if not os.path.isfile(mut_file_name):
            return f"Error, no mut file for variant {variant_path}"
        else:
            command_mut = f"pyDock3 {mut_file_name} opra"
            sp.run(command_mut, shell=True)

    return f"Successfully created oda file for variant {variant_path}"