import pandas as pd
from data_retrievel_and_feature_extraction import uniprot_info as uni


def add_aa_position_to_df_according_to_seq(df: pd.DataFrame, seq: str) -> pd.DataFrame:
    subseq = seq
    positions_list = list(range(1, len(seq)+1)) #This list is consists of positions corresponding with the aa positions in the sequence.
    output_df = pd.DataFrame()
    # create subdf, which is df with "aa_pos" column, initiated to null.
    subdf = df.copy()
    subdf['aa_pos'] = None

    for aa, pos in zip(seq, positions_list):
        if not is_next_aa_the_same(subseq):
            df_for_aa_with_position = get_df_rows_for_aa_with_position_by_aa_seq(subdf, pos, aa)
            # add intermidiate_df to output_df
            output_df = pd.concat([output_df, df_for_aa_with_position])
            subdf = subdf[~subdf["aa_pos"].notna()]  # remove rows that correpond with df_for_aa_with_position
            subseq = subseq[1:]
        elif is_next_aa_the_same(subseq):
            # I will only add by the DNA, and will not treat the next aa that is the same! It will loop again for this.
            df_for_aa_with_position = get_df_rows_for_aa_with_position_with_multiple_concequitive(subdf, pos) # Intermidiate_df will have only the rows with position.
            output_df = pd.concat([output_df, df_for_aa_with_position])
            subdf = subdf[~subdf["pos_hg19"].isin(df_for_aa_with_position["pos_hg19"])] # remove rows that correpond with df_for_aa_with_position
            subseq = subseq[1:]
    df = output_df
    return df


def get_df_rows_for_aa_with_position_with_multiple_concequitive(subdf: pd.DataFrame, pos) -> pd.DataFrame:
    """
    This method adds the pos to the aa by the addition in the dna column.
    This will only add pos in position column for the aa in this pos.
    Returns:
        Dataframe like input df, with only rows
    """
    first_dna_value = None
    for index, row in subdf.iterrows():
        current_dna_value = row['pos_hg19'] # get dna value.
        if first_dna_value is None:
            first_dna_value = row['pos_hg19']
            subdf.at[index, "aa_pos"] = pos
        elif current_dna_value - first_dna_value >= 3: # It's a different aa.
            break
        else:
            subdf.at[index, "aa_pos"] = pos # add position to "aa_pos" column in this row.
    subdf = subdf[subdf['aa_pos'].notna()]
    return subdf


def get_df_rows_for_aa_with_position_by_aa_seq(subdf: pd.DataFrame, pos, aa):
    print(aa)
    for index, row in subdf.iterrows():
        if row["aa_wt"] == aa:
            subdf.at[index, "aa_pos"] = pos
        else:
            break
    subdf = subdf[subdf['aa_pos'].notna()]
    return subdf


def is_next_aa_the_same(seq: str) -> bool:  # tested and validated
    if len(seq) <= 1:
        return False
    elif seq[0] == seq[1]:
        return True
    return False


def get_consequitive_num(subseq: str) -> int:  # Not used
    aa = subseq[0]
    count = 0
    for i in subseq:
        if subseq[i] != aa:
            break
        else:
            count+=1
    return count


def test_fit_to_seq(df, seq):
    problematic_positions = []
    for pos, letter in enumerate(seq, 1):
        # Go over the position and aa_wt in every row and compare to pos and letter.
        matching_rows = df.loc[df["aa_pos"] == pos]

        # Check if the 'aa_wt' column in the matching rows contains 'letter'
        for index, row in matching_rows.iterrows():
            if row['aa_wt'] == letter:
                continue
            else:
                print(f"Error {index}: aa_wt does not match the sequence at position {pos}.")
                problematic_positions.append(pos)
    return set(problematic_positions)


def test_all_positions_in_df(df, seq):
    all_positions_in_sequence = set(range(1, len(seq) + 1))
    positions_in_dataframe = set(df['aa_pos'])
    # Find the positions in the sequence that are not in the DataFrame
    positions_not_in_dataframe = all_positions_in_sequence - positions_in_dataframe

    return positions_not_in_dataframe


def add_aa_position_to_df_reverse(df: pd.DataFrame) -> pd.DataFrame:
    """This function adds the amino acid position to the dataframe, under new column 'aa_pos'
    Method:
    Add a new column 'aa_pos' to the dataframe
    Iterate through the "pos_hg19" column, and for each increasment of 3, add the value to the new column
    """
        # Initialize variables to keep track of the last value and the incremental number
    last_value = None
    incremental_number = 0
    number_column = []

    # Iterate through the pos_hg19 column and calculate the numbers column
    for value in df['pos_hg19']:
        if last_value is None:
            last_value = value
            incremental_number = 1
        elif last_value - value >= 3:
            incremental_number += 1
            last_value = value
        else:
            pass
        number_column.append(incremental_number)

    # Add the numbers column to the DataFrame
    df['aa_pos'] = number_column
    return df


def encode_amino_acids(df):
    df = df.sort_values(by='aa_pos').reset_index(drop=True)  # Sort the DataFrame by 'aa_pos' and reset the index
    aa_seq = ''

    for index, row in df.iterrows():
        aa_wt = row['aa_wt']
        aa_seq += aa_wt

    return aa_seq


if __name__ == "__main__":
    reverse = False
    gene = 'COL4A3'
    if reverse:
        revel_data = pd.read_csv(f"{gene}_reverse.csv", header=None)
    else:
        revel_data = pd.read_csv(f"{gene}.csv", header=None)
    revel_data.columns = ['chr', 'pos_hg19', 'pos_grch38', 'ref_na', 'alt_na', 'aa_wt', 'aa_mut', 'revel_score',
                          "transcript_id"]
    seq = uni.get_sequence(gene)

    if reverse:
        revel_data = revel_data[::-1].reset_index(drop=True)
        revel_data = add_aa_position_to_df_reverse(revel_data)
    else:
        revel_data = add_aa_position_to_df_according_to_seq(revel_data, seq)

