# This function creates the mutated protein usind FoldX
# $! is the pdb file.
# $2 is the mutation position.
# $3 is the outpur directory

input_string=$2
modified_string="${input_string:0:1}A${input_string:1}"

# Create the output directory if it doesn't exist
mkdir -p "$3"

# Run the FoldX command with the output path
foldx5 --command=PositionScan --pdb=$1 --positions=$modified_string --out-path="$3"
