import pandas as pd
import os

def get_revel_file(chromosome, position):
    """This function gets the chromosome and position of a variant, and returns the right revel file to get the revel score from.
    The revel file names are in the format: "revel_grch38_chrom_{num}_{first_position}_{last_position}.csv"
    """
    files = os.listdir("C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\benchmarking\\REVEL\\revel_prediction_files")
    for file in files:
        if file.split("_")[3] == str(chromosome) and int(file.split("_")[4]) <= position <= int(file.split("_")[5].split(".")[0]):
            return file
    print("No revel file found for chromosome " + str(chromosome) + " and position " + str(position) + ".")
    return None


def get_revel_score(file_path, pos, ref, alt, aaref, aaalt) -> int:
    """This function receives a path to a revel file, and the position, ref, alt, aaref and aaalt of a variant, and returns the revel score of the variant. If the variant is not in the file, it returns None.
    """
    revel_df = pd.read_csv(file_path, low_memory=False)
    revel_score = revel_df[(revel_df["hg19_pos"] == pos) & (revel_df["ref"] == ref) & (revel_df["alt"] == alt) & (revel_df["aaref"] == aaref) & (revel_df["aaalt"] == aaalt)]["REVEL"]
    if len(revel_score) > 0:
        return revel_score.values[0]
    else:
        return None


if __name__ == "__main__":
    # For each dvd file, get the revel score for each variant in the dataset.
    # The revel scores are in the "C:\Users\InbarBlech\PycharmProjects\Thesis\benchmarking\REVEL\revel_prediction_files" directory.

    files = os.listdir("C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\Data\\all_dvd_data_with_genomic_info\\")

    revel_scores = pd.DataFrame(columns=["variant", "gene", "revel_score"])

    # initiate a dataframe that will contain the errors and successes of the process, and the variables that were used.
    errors_df = pd.DataFrame(
        columns=["file", "variant", "gene", "error", "gene_variant_wt", "gene_variant_mut", "protein_variant_wt",
                 "protein_variant_mut", "chromosome", "position"])

    for file in files:
        # Read the dvd file
        df = pd.read_csv(
            "C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\Data\\all_dvd_data_with_genomic_info\\" + file)
        # Get the revel score for each variant in the dataset from the right revel file in the revel_prediction_files directory, and add it to the revel_scores dataframe.
        df['chromosome'] = df["variation"].apply(lambda x: x.split(":")[0] if ":" in x else None)
        for index, row in df.iterrows():
            chromosome = row["chromosome"]
            if chromosome != "X" and chromosome != "Y":
                if int(chromosome) < 10:
                    chromosome = "0" + chromosome
            position = row["pos"]
            protein_variant_wt = row["variant"][0]
            protein_variant_mut = row["variant"][-1]
            gene_variant_wt = row["variation"][-3]
            gene_variant_mut = row["variation"][-1]
            # Get the right revel file according to the chromosome and position of the variant.
            revel_file = get_revel_file(chromosome, position)
            # Extract the revel score from the revel file, according to the position and protein variant of the variant.
            if revel_file is not None:
                revel_score = get_revel_score(
                    "C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\benchmarking\\REVEL\\revel_prediction_files\\" + revel_file,
                    position, gene_variant_wt, gene_variant_mut, protein_variant_wt, protein_variant_mut)
                if revel_score is not None:
                    # Add the revel score to the revel_scores dataframe, without using append method
                    new_row = {"variant": row["variant"], "gene": row["gene"], "revel_score": revel_score}
                    revel_scores = pd.concat([revel_scores, pd.DataFrame([new_row])], ignore_index=True)
                    print("Revel score for variant " + row["variant"] + " in file " + revel_file + " is " + str(
                        revel_score) + ".")
                    # Add the success to the errors dataframe
                    new_row = {"file": file, "variant": row["variant"], "gene": row["variation"], "error": "Success",
                               "gene_variant_wt": gene_variant_wt, "gene_variant_mut": gene_variant_mut,
                               "protein_variant_wt": protein_variant_wt, "protein_variant_mut": protein_variant_mut,
                               "chromosome": chromosome, "position": position}
                    errors_df = pd.concat([errors_df, pd.DataFrame([new_row])], ignore_index=True)
                else:
                    print("No revel score found for variant " + row["variant"] + " in file " + revel_file + ".")
                    # Add the error to the errors dataframe
                    new_row = {"file": file, "variant": row["variant"], "gene": row["variation"],
                               "error": "No revel score found", "gene_variant_wt": gene_variant_wt,
                               "gene_variant_mut": gene_variant_mut, "protein_variant_wt": protein_variant_wt,
                               "protein_variant_mut": protein_variant_mut, "chromosome": chromosome,
                               "position": position}
                    errors_df = pd.concat([errors_df, pd.DataFrame([new_row])], ignore_index=True)
            else:
                print("No revel file found for chromosome " + str(chromosome) + " and position " + str(position) + ".")
                # Add the error to the errors dataframe
                new_row = {"file": file, "variant": row["variant"], "gene": row["variation"],
                           "error": "No revel file found", "gene_variant_wt": gene_variant_wt,
                           "gene_variant_mut": gene_variant_mut, "protein_variant_wt": protein_variant_wt,
                           "protein_variant_mut": protein_variant_mut, "chromosome": chromosome, "position": position}
                errors_df = pd.concat([errors_df, pd.DataFrame([new_row])], ignore_index=True)
        print("Finished file " + file + ".")
        # Save the revel scores dataframe to a file.
        revel_scores.to_csv(
            "C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\benchmarking\\REVEL\\revel_scores_for_all_dvd.csv",
            index=False)
        # Save the errors dataframe to a file.
        errors_df.to_csv(
            "C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\benchmarking\\REVEL\\errors_for_revel_scores_for_all_dvd.csv",
            index=False)
