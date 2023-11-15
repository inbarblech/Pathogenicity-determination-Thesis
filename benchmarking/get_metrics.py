import csv
import os
import random
import shutil

import pandas as pd


def change_file_names():
    # Change the name of all csv files in the folder
    folder_path = "/home/inbar/mutpred/"

    # List all the CSV files in the folder
    csv_files = [f for f in os.listdir(folder_path) if f.endswith('.csv')]

    # Shuffle the list of files to randomize their order
    random.shuffle(csv_files)

    # Iterate through the shuffled files and rename them to numbers
    for i, csv_file in enumerate(csv_files, start=1):
        new_name = f"{i}.csv"
        old_path = os.path.join(folder_path, csv_file)
        new_path = os.path.join(folder_path, new_name)
        os.rename(old_path, new_path)
        print(f'Renamed {csv_file} to {new_name}')

    print("File renaming complete.")


def add_gene_name_column_to_csv(file_path):
    """This function adds a gene name column to the csv file in the given folder.
    The gene name is taken from the file name (the first part of the file name)."""
    # Get the gene name from the file name
    gene_name = file_path.split("\\")[-1].split('_')[0]
    print(gene_name)

    # Create a new CSV file
    with open(file_path, 'r') as infile:
        # convert to dataframe
        df = pd.read_csv(infile)
        # add gene name column
        df['gene'] = gene_name
        # convert back to csv
        df.to_csv(file_path, index=False)


def create_one_csv_from_all_csv_files():
    # Create one CSV file from all the CSV files in the folder
    # folder_path = "/home/inbar/predictions/EVE/"
    folder_path = "C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\benchmarking\\REVEL\\"

    # List all the CSV files in the folder
    csv_files = [f for f in os.listdir(folder_path) if f.endswith('.csv')]

    # Create a new CSV file
    with open(f"{folder_path}all_REVEL_predictions.csv", 'w') as outfile:
        for i, csv_file in enumerate(csv_files, start=1):
            with open(os.path.join(folder_path, csv_file)) as infile:
                for line in infile:
                    outfile.write(line)
            print(f"Added {csv_file} to all_REVEL_predictions.csv")


if __name__ == "__main__":
    genes = {"MYO7A", "FGFR1", "WFS1", "COL2A1", "COL4A3", "COL4A5"}
    for gene in genes:
        add_gene_name_column_to_csv(f"C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\benchmarking\\REVEL\\{gene}_revel_with_pos.csv")
    create_one_csv_from_all_csv_files()