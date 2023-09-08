# Writing to file
with open("C:\\Users\\InbarBlech\\Downloads\\uniprots_without_msa.txt", 'w') as f:
    for gene in uniprot_list_without_msa:
        f.write(f"{gene}\n")


# saving figure
