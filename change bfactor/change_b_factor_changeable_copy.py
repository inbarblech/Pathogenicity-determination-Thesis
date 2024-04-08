import pandas as pd
import os

"""
This script changes the b-factor column in the pdb file with the b-factor values in the csv file
The csv file contains one column with the b-factor values, and the pdb file contains the atom rows alone.
The values in the csv file are 00.00 or 10.00 or 20.00 or 30.00 and so on.

Before run:
-   Change the directory_path and pdb_files_path to the relevant paths, or change the folder variable to the relevant folder.
-   First run a function that creates the b_factor series, and then run this function.
    For example for such function, look at "structural_analysis_of_predictor_preformance_using_bfactor_change_file_setup.py"

Make sure you have created the directories:
    - "b_factor_pdbs" in the directory_path (Inside the directory_path)
    - There are edited pdb files in the pdb_files_path, with only the atom rows.
    - the b_factor and pdb_files directories inside the folder in the predictions_vs_real directory.

This function was written to change the b-factor according to the predictions of the model, and it's input is either "pathogenic" or "benign",
for the distinction between TP, TN, FP and FN.
If you do not want to use this distinction, change the b_factor_purpose to "all" and remove the if statements in the for loop.
It was not tested for this purpose, so make sure to check the output.
"""

genes = ["COL2A1", "COL4A3", "COL4A5", "GJB2", "FGFR1", "MYO7A", "SLC26A4", "WFS1"]
folder = "gene_specific_predictors"
b_factor_purpose = ["benign", "pathogenic"]

for gene in genes:
    directory_path = f"C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\predictions_vs_real\\{folder}\\{gene}\\performance_using_b_factor\\"
    pdb_files_path = f"C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\predictions_vs_real\\{folder}\\{gene}.pdb"
    # For every file in the directory, change the b-factor column in the pdb file
    # with the b-factor values in the csv file
    items = os.listdir(directory_path)
    files = [file for file in items if file.endswith('.csv')]
    for purpose in b_factor_purpose:
        for file in files:
            if purpose == "pathogenic":
                # Only use files with the word pathogenic in them
                if "pathogenic" not in file:
                    continue
            if purpose == "benign":
                # Only use files with the word benign in them
                if "benign" not in file:
                    continue
            gene = file.split('.')[0].split('_')[0]
            pdb_path = f'{pdb_files_path}'
            csv_b_factor_file_path = f'{directory_path}{file}'
            pdb_output_path = f'C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\predictions_vs_real\\{folder}\\performance_in_structure_b_factor\\{gene}{purpose}.pdb'
            df = pd.read_csv(csv_b_factor_file_path)

            with open(pdb_path, 'r') as f:
                lines_f = f.readlines()
                new_lines = []
                for index, line in enumerate(lines_f):
                    new_line = line[:61] + str(df.loc[index, 'b_factor']) + line[66:]
                    new_lines.append(new_line)

            with open(pdb_output_path, 'w') as f:
                f.writelines(new_lines)