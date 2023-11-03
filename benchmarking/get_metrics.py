import os
import random
import shutil

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


def create_one_csv_from_all_csv_files():
    # Create one CSV file from all the CSV files in the folder
    folder_path = "/home/inbar/predictions/VEST4/"

    # List all the CSV files in the folder
    csv_files = [f for f in os.listdir(folder_path) if f.endswith('.csv')]

    # Create a new CSV file
    with open(f"{folder_path}all.csv", 'w') as outfile:
        for i, csv_file in enumerate(csv_files, start=1):
            with open(os.path.join(folder_path, csv_file)) as infile:
                for line in infile:
                    outfile.write(line)
            print(f"Added {csv_file} to all_vest4_predictions.csv")


if __name__ == "__main__":
    create_one_csv_from_all_csv_files()