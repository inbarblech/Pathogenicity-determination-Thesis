import subprocess as sp
import os


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


if __name__ == '__main__':
    clean_file_names('/home/inbar/DVDdata/Benign/')
    clean_file_names('/home/inbar/DVDdata/Pathogenic/')
