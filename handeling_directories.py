import shutil
import subprocess as sp
import os


def check_size_of_directory(path_to_gene_folder: str):
    """Checks the amount of files in each folder in the given path, including directories.
    If the amount is less than 26, return the folder name, and the amount of files in it.
    Args:
        path_to_gene_folder (str): The path to the gene folder.
    return:
        errors (list): A list of the folders that have less than 26 files in them.
    """
    errors = []

    for folder in os.listdir(path_to_gene_folder):
        folder_path = os.path.join(path_to_gene_folder, folder)

        if os.path.isdir(folder_path):  # Check if the path is a directory
            if 26 < len(os.listdir(folder_path)) < 0:
                errors.append(f"{folder} ({len(os.listdir(folder_path))} \n")

    return errors


def run_check_size_on_gene_folders(path):
    """Runs a function on all gene folders in path"""
    errors = []
    os.chdir(path)
    items_in_folder = os.listdir(path)
    gene_folders_paths = [os.path.join(path, item) for item in items_in_folder if
                          os.path.isdir(os.path.join(path, item))]
    for gene_folder in gene_folders_paths:
        print("Checking size of directory: ", gene_folder)
        errors.append(check_size_of_directory(gene_folder))

    with open(f"{path}/directories_without_all_files.txt", "w") as f:
        for error in errors:
            f.write(f"{error}\n")


def check_if_all_folders_contain_a_type_of_file(path_to_gene_folder: str, file_wt, file_mut,
                                                output_file_name_wt, output_file_name_mut):
    """Checks if all the folders in the given path contain the give file.
    If not, writes the folder name to a file.
    Works for wt and mut separately.
    For example, if we want to see if all the folders contain a file called "wt.txt", we will call the function with
    file_wt = "wt.pdb" and file_mut = "mut.pdb".
    The function will write to a file called "errors_wt.txt" all the folders that do not contain a file called "wt.txt".
    Args:
        path_to_gene_folder (str): The path to the gene folder.
        file_wt (str): The name of the file to check if it exists in all the folders.
        file_mut (str): The name of the file to check if it exists in all the folders.
        output_file_name_wt (str): The name of the file to write the errors for wt.
        output_file_name_mut (str): The name of the file to write the errors for mut.
    """
    errors_wt = []
    errors_mut = []

    for folder in os.listdir(path_to_gene_folder):
        folder_path = os.path.join(path_to_gene_folder, folder)

        if os.path.isdir(folder_path):  # Check if the path is a directory
            if not os.path.isfile(os.path.join(folder_path, file_wt)):
                errors_wt.append(folder)
            if not os.path.isfile(os.path.join(folder_path, file_mut)):
                errors_mut.append(folder)

            if len(os.listdir(folder_path)) == 0:
                # Remove empty folders from the lists
                if folder in errors_wt:
                    errors_wt.remove(f"folder does not contain file for wt: {folder}")
                if folder in errors_mut:
                    errors_mut.remove(f"folder does not contain file for mut: {folder}")

    with open(f"{path_to_gene_folder}/{output_file_name_wt}", "w") as f:
        for error in errors_wt:
            f.write(f"{error}\n")

    with open(f"{path_to_gene_folder}/{output_file_name_mut}", "w") as f:
        for error in errors_mut:
            f.write(f"{error}\n")


def delete_folder(path):
    """Deletes the folder in the given path, even if it's not empty."""
    sp.run(f"rm -rf {path}", shell=True)


def clean_file_names(directory):
    for root, dirs, files in os.walk(directory):
        for file_name in files:
            original_path = os.path.join(root, file_name)
            parts = file_name.split('.')
            if len(parts) > 4:
                new_file_name = '.'.join(parts[:4]) + '.csv'
                new_path = os.path.join(root, new_file_name)
                os.rename(original_path, new_path)
                print(f'Renamed "{file_name}" to "{new_file_name}"')


def change_file_name(folder_path: str, to_replace: str, replace_with: str):
    """Changes the file name of the given file path to the given new name.
    Args:
        folder_path (str): The path to the folder containing the file.
        to_replace (str): The string to replace in the file name.
        replace_with (str): The string to replace with in the file name.
    """

    for filename in os.listdir(folder_path):
        if filename.endswith(to_replace):
            new_filename = filename.replace(to_replace, replace_with)
            old_path = os.path.join(folder_path, filename)
            new_path = os.path.join(folder_path, new_filename)
            os.rename(old_path, new_path)
            print(f'Renamed: {filename} -> {new_filename}')


def run_on_gene_folders(path, list_of_genes_to_run_on=None):
    """Runs a function on all gene folders in path"""
    os.chdir(path)
    # get list of gene folders paths
    if list_of_genes_to_run_on is None:
        items_in_folder = os.listdir(path)
    else:
        items_in_folder = list_of_genes_to_run_on
    gene_folders_paths = [os.path.join(path, item) for item in items_in_folder if
                          os.path.isdir(os.path.join(path, item))]
    for gene_folder in gene_folders_paths:
        run_change_name_on_folder(gene_folder)


def count_number_of_variant_folders(path_to_folder):
    """Counts the number of folders and sub-folders in the given folder."""
    counter = 0
    for root, dirs, files in os.walk(path_to_folder):
        for dir in dirs:
            counter += 1
    return counter


def run_change_name_on_folder(path_to_gene_folder):
    """Runs on all variant folders in path_to_folder
    Usage:
        run_change_name_on_folder(path_to_gene_folder)
    """
    errors = []
    # change directory to the folder containing all the gene folders
    os.chdir(path_to_gene_folder)
    # get list of variant folders paths
    all_items_in_gene_folder = os.listdir(path_to_gene_folder)

    # get list of variant folder paths
    variant_folders_paths = [os.path.join(path_to_gene_folder, item) for item in all_items_in_gene_folder if
                             os.path.isdir(os.path.join(path_to_gene_folder, item))]

    print(f"variant_folders_paths: {variant_folders_paths}")
    # create mutation to all variants in all genes, using the function add_mut.create_mut_using_foldx
    for variant_folder in variant_folders_paths:
        change_file_name(variant_folder, 'pdb.oda', 'oda.pdb')
        print(f"variant_folder: {variant_folder}")


def check_two_file_type_exists_in_variant_path(variant_path: str, type: str) -> bool:
    """Checks if the oda file exists for the given variant.
    Args:
        variant_path (str): The path to the variant folder.
    Returns:
        bool: True if the oda file exists, False otherwise.
    """
    # change the working directory to the variant folder
    os.chdir(variant_path)
    # list all files in the variant folder that end with the given type
    files = [file for file in os.listdir(variant_path) if file.endswith(type)]
    # if there is no file that ends with the given type, return False
    if len(files) == 0 or len(files) == 1:
        return False
    elif len(files) == 2:
        return True
    else:
        raise False


def check_file_type_exists_in_variant_path(variant_path: str, type: str) -> bool:
    """Checks if the oda file exists for the given variant.
    Args:
        variant_path (str): The path to the variant folder.
    Returns:
        bool: True if the oda file exists, False otherwise.
    """
    # change the working directory to the variant folder
    os.chdir(variant_path)
    # list all files in the variant folder that end with the given type
    files = [file for file in os.listdir(variant_path) if file.endswith(type)]
    # if there is no file that ends with the given type, return False
    if len(files) == 0:
        return False
    else:
        return True


def copy_file_from_one_directory_to_all_folders_in_path(file_path: str, path_of_directories: str):
    """Copies the given file to all folders in the given path.
    Args:
        file_path (str): The path to the file to copy.
        path_of_directories (str): The path to the folders to copy the file to. (Gene folder)
    """
    source_directory = os.path.dirname(file_path)
    for root, dirs, files in os.walk(path_of_directories):
        for directory in dirs:
            # If the directory is the same as the file path, skip it
            current_directory = os.path.join(root, directory)
            # If the current directory is the same as the source directory, skip it
            if current_directory == source_directory:
                continue
            shutil.copy(file_path, os.path.join(path_of_directories, directory))


if __name__ == "__main__":


    files_with_no_oda = []
    files_with_no_sasa = []

    # # copy wt oda file to all variant folders in genes
    # paths_of_directories = ["/home/inbar/variants/Benign_for_gene_specific/SLC26A4/",
    #                         "/home/inbar/variants/Benign_for_gene_specific/WFS1/",
    #                         "/home/inbar/variants/Benign_for_gene_specific/COL4A5/",
    #                         "/home/inbar/variants/Benign_for_gene_specific/COL4A3/",
    #                         "/home/inbar/variants/Benign_for_gene_specific/COL2A1/",
    #                         "/home/inbar/variants/Benign_for_gene_specific/MYO7A/",
    #                         "/home/inbar/variants/Benign_for_gene_specific/FGFR1/"]
    # file_paths = ["/home/inbar/variants/Benign_for_gene_specific/SLC26A4/SLC26A4_O43511_G5R/AF_O43511.pdb.oda",
    #                    "/home/inbar/variants/Benign_for_gene_specific/WFS1/WFS1_O76024_P7A/AF_O76024.pdb.oda",
    #                    "/home/inbar/variants/Benign_for_gene_specific/COL4A5/COL4A5_P29400_I43V/AF_P29400.pdb.oda",
    #                     "/home/inbar/variants/Benign_for_gene_specific/COL4A3/COL4A3_Q01955_Q75T/AF_Q01955.pdb.oda",
    #                    "/home/inbar/variants/Benign_for_gene_specific/COL2A1/COL2A1_P02458_Q8P/AF_P02458.pdb.oda",
    #                    "/home/inbar/variants/Benign_for_gene_specific/MYO7A/MYO7A_Q13402_R388K/AF_Q13402.pdb.oda",
    #                    "/home/inbar/variants/Benign_for_gene_specific/FGFR1/FGFR1_P11362_M1T/AF_P11362.pdb.oda"]
    # for file_path, paths_of_directory in zip(file_paths, paths_of_directories):
    #     copy_file_from_one_directory_to_all_folders_in_path(file_path, paths_of_directory)
    #
    # file_paths = ["/home/inbar/variants/Benign_for_gene_specific/SLC26A4/SLC26A4_O43511_G5R/AF_O43511.asa.pdb",
    #               "/home/inbar/variants/Benign_for_gene_specific/WFS1/WFS1_O76024_P7A/AF_O76024.asa.pdb",
    #               "/home/inbar/variants/Benign_for_gene_specific/COL4A5/COL4A5_P29400_I43V/AF_P29400.asa.pdb",
    #               "/home/inbar/variants/Benign_for_gene_specific/COL4A3/COL4A3_Q01955_Q75T/AF_Q01955.asa.pdb",
    #               "/home/inbar/variants/Benign_for_gene_specific/COL2A1/COL2A1_P02458_Q8P/AF_P02458.asa.pdb",
    #               "/home/inbar/variants/Benign_for_gene_specific/MYO7A/MYO7A_Q13402_R388K/AF_Q13402.asa.pdb",
    #               "/home/inbar/variants/Benign_for_gene_specific/FGFR1/FGFR1_P11362_M1T/AF_P11362.asa.pdb"]
    # for file_path, paths_of_directory in zip(file_paths, paths_of_directories):
    #     copy_file_from_one_directory_to_all_folders_in_path(file_path, paths_of_directory)

    # # Check if all variant folders in genes have 2  file
    # for gene in ["COL2A1", "COL4A3", "COL4A5", "FGFR1", "MYO7A", "SLC26A4", "WFS1"]:
    #     gene_path = f"/home/inbar/variants/Benign_for_gene_specific/{gene}/"
    #     variant_paths = [os.path.join(gene_path, item) for item in os.listdir(gene_path) if
    #                      os.path.isdir(os.path.join(gene_path, item))]
    #     for path in variant_paths:
    #         if check_two_file_type_exists_in_variant_path(path, "pdb.oda"):
    #             pass
    #             # Change the name of the file from pdb.oda to oda.pdb
    #             # change_file_name(path, 'pdb.oda', 'oda.pdb')
    #         else:
    #             files_with_no_oda.append(path)
    #         if check_two_file_type_exists_in_variant_path(path, "asa.pdb"):
    #             pass
    #         else:
    #             files_with_no_sasa.append(path)

    # with open("C:\\Users\\InbarBlech\\Downloads\\asamissingfiles210923.txt", 'w') as f:
    #     for gene in files_with_no_sasa:
    #         f.write(f"{gene}\n")
    # with open("C:\\Users\\InbarBlech\\Downloads\\odamissingfiles210923.txt", 'w') as f:
    #     for gene in files_with_no_oda:
    #         f.write(f"{gene}\n")
    # print(f"files_with_no_sasa: {files_with_no_sasa}")
    # # print(f"files_with_no_oda: {files_with_no_oda}")
    # print(f"len(files_with_no_sasa): {len(files_with_no_sasa)}")
    # # print(f"len(files_with_no_oda): {len(files_with_no_oda)}")
    # #
    # # Check overlap between files_with_no_sasa and files_with_no_oda
    # overlap = set(files_with_no_sasa).intersection(files_with_no_oda)
    # # print(f"overlap: {overlap}")
    # print(f"len(overlap): {len(overlap)}")
    #
    # # print the files that are in files_with_no_sasa but not in files_with_no_oda
    # files_with_no_oda_but_with_sasa = set(files_with_no_sasa).difference(files_with_no_oda)
    # print(f"files_with_no_oda_but_with_sasa: {files_with_no_oda_but_with_sasa}")
    #
    # for gene in ["COL2A1", "COL4A3", "COL4A5", "FGFR1", "MYO7A", "SLC26A4", "WFS1"]:
    #     gene_path = f"/home/inbar/variants/Benign_for_gene_specific/{gene}/"
    #     variant_paths = [os.path.join(gene_path, item) for item in os.listdir(gene_path) if
    #                      os.path.isdir(os.path.join(gene_path, item))]
    #     for path in variant_paths:
    #         # Change the name of the file from pdb.oda to oda.pdb
    #         # change_file_name(path, 'pdb.oda', 'oda.pdb')
    #         # Check if there is any file that ends with pdb.oda
    #         if not check_two_file_type_exists_in_variant_path(path, "oda.pdb"):
    #             print(f"pdb.oda file exists in {path}")

    run_change_name_on_folder("/home/inbar/variants/Benign_for_gene_specific/GJB2/")





