import pandas as pd
## Write function that creates a txt file in this format:
# VAR1 np_.. A24P
# VAR2 np_.. A2545M
# VAR3 np_.. A453T
# For each variant in the input csv file, create a line in the txt file.


def create_txt_file_from_csv(csv_file_path: str, crv_file_path: str, np: str):
    """Create a txt file from the given csv file.
    Each line in the txt file will be in this format:
    VAR1 NP_.. A24P
    VAR2 NP_.. A2545M
    VAR3 NP_.. A453T
    """
    # Read the csv file
    df = pd.read_csv(csv_file_path)

    for index, row in df.iterrows():
        # append the data for this row in the format: VAR{index} np_.. {variant}
        with open(crv_file_path, 'a') as file:
            file.write(f'TR{index} {np} {row["variant"]}\n')

    print(f'Created txt file in {crv_file_path}')


if __name__ == "__main__":
    genes = ["GJB2"]
    file_with_variants_name = None  # Initiated inside the for loop
    gene_np_ids_dict = {"SLC26A4": "NP_000432.1", "COL4A3": "NP_000082.2", "COL4A5": "NP_203699.1", "COL2A1": "NP_001835.3", "WFS1": "NM_006005.3", "FGFR1": "NP_075598.2", "MYO7A": "NP_000251.3", "GJB2": "NP_003995"}
    for gene in genes:
        file_with_variants_name = f"{gene}_predictions_LOPO_XGB.csv"
        create_txt_file_from_csv(f"C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\predictions_vs_real\\{gene}\\{file_with_variants_name}", f"C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\benchmarking\\VEST4\\VEST4_input_{gene}.crv", f"{gene_np_ids_dict[gene]}")
