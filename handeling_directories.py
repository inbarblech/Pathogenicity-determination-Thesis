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
    """Runs on all variant folders in path_to_folder"""
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



if __name__ == "__main__":
    run_on_gene_folders("/home/inbar/variants/Pathogenic/")
    run_on_gene_folders("/home/inbar/variants/Benign/")