# This function produces the files of the different features using the various function.

# VARIABLES:
# $wild_type_pdb is the name of file (or path if in different directory) of wt protein pdb.
# $mutated_pdb is the name of file (or path if in different directory) of mutated protein pdb.

# NOTE TO SELF : this function should be run from a directory that is specific to the residue/mutation, since the files will be saved without any recognition (In the name there is just wt or mut)

# Optimal docking area
pyDock3 $wild_type_pdb oda > wt.ODAtab
pyDock3 $mutated_pdb oda > mut.ODAtab
grep '.$mut_loc'*.ODAtab
oda_mut=
oda_wt=
#Here need to compare between the values of wild_type and mutated_protein and add them to variables
delta_oda=$oda_mut-$oda_wt #need to create those variables

# Opra
pyDock3 $wild_type_pdb opra > wt.opra.pdb
pyDock3 $mutated_pdb opra > mut.opra.pdb
awk '{ if ($11 <= -1.0) print ;)' $opra_file
opra_wt=
opra_mut=
#Here need to compare between the values of wild_type and mutated_protein and add them to variables
delta_opra=$opra_mut-$opra_wt #need to creat these variables

# Solvant accesible surface area (SASA)
dr_sasa -m 0 -i $wild_type_pdb > wt.asa.pdb
dr_sasa -m 0 -i $mutated_pdb > mut.asa.pdb
grep 'ALA A  19' asa.pdb | awk '{ SUM += $11} END { print SUM }'
sasa_wt=
sasa_mut=
#Here need to compare between the values of wild_type and mutated_protein and add them to variables
delta_sasa=$sasa_mut-$sasa_wt #need to creat these variables

# Here need to create one line of all variables and add it to an existing csv file of all features.
