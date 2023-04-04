# This script adds a variant to the csv, including all the features
# Input: uniprot code (Ex. P51841), gene name (Ex. GUCY2F), variant (Ex. G690S)
# Output: CSV file with added variant features.

# Requirments: All necessary programs: Dr_SASA, pyDoch, FoldX5, access to internet for fetching PDB files, consurf, p2rank.

# making a new directory for the mutation
# $1 is uniprot, $2 is Gene name, $3 is mutation
cd
cd variants 
mkdir "$1_$2"
cd "$1_$2"
mkdir "$3" # This is the mutation
cd "$3"

# Fetching the PDB file
wget http://alphafold.ebi.ac.uk/files/AF-"$1"-F1-model_v4.pdb


