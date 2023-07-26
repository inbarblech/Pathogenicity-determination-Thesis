import os
import uniprot_info as uni
import subprocess as sp
import dvd_data as dvd
import handeling_directories as hd


def check_aa_in_sequence(aa, seq, pos):
    """Returns True if the given amino acid is in the given sequence, in the right position.
    """
    if aa == seq[pos-1]:
        return True
    else:
        return False


def create_path_to_file(path, file_name):
    """Creates a file name at the given path.
    if the path doesn't exist, creates it.
    if the path exists, returns False."""
    new_path = os.path.join(path, file_name)
    if not os.path.exists(new_path):
        os.makedirs(new_path)
        print(f"""-----------------
        path {file_name} created successfully.""")
        return new_path
    else: # if the path exists
        error_message = f"Path {new_path} already exists."
        return "Error"


def make_file(path, file_name):
    """Makes an empty file at the given path."""
    with open(path, 'w') as f:
        f.write('')


def get_structure_af(gene_id):
    """Returns the PDB file for the given mutation."""
    # Save the PDB file to disk
    sp.run(f"wget http://alphafold.ebi.ac.uk/files/AF-{gene_id}-F1-model_v4.pdb", shell=True)
    # Return the path to the downloaded PDB file
    return f"{gene_id}.pdb"


def create_mutation_using_foldx(mut, gene_id):
    """creates a mutation using foldx"""
    mut_with_a = f"{mut[0]}A{mut[1:]}"
    command = f"foldx --command=PositionScan --pdb=AF_{gene_id}.pdb --positions={mut_with_a}"
    sp.run(command, shell=True)


if __name__ == "__main__":
    path = "/home/inbar/variants/Benign/POLR1B/POLR1B_Q9H9Y6_N9S/"
    sp.run()

