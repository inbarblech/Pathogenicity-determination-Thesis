import pandas as pd
from data_retrievel_and_feature_extraction import uniprot_info as uni



# def add_aa_position_to_df_according_to_seq(df: pd.DataFrame, seq: str) -> pd.DataFrame:
#     """This function adds the amino acid position to the dataframe, under new column 'aa_pos'
#     Method:
#     Add a new column 'aa_pos' to the dataframe
#     Iterate through the "pos_hg19" column, and for each increasment of 3, add the value to the new column
#     This version of the function also checks that the amino acid chain is corresponding with the input sequence.
#     """
#     # Initialize variables to keep track of the last value and the incremental number
#     last_value = None
#     incremental_number = 0
#     number_column = []
#     aa_list = [aa for aa in seq]
#     aa = None
#     consecutives = 0
#     sub_dataframe = None
#
#     for aa in range(len(seq) - 1):
#         incremental_number =+ 1
#         while is_next_aa_the_same():
#             # do something
#             # Add to consecutives
#             consecutives +=1
#         if consecutives !=0:
#             # Treat the counting accordingly - There are more than 1 consecutives aa that's the same.
#             # Then do consecutives = 0
#             consecutives = 0
#         else: # consecutives is zero, meaning you only need to treat one letter, add one to the counting.
#
#     # Add the numbers column to the DataFrame
#     df['aa_pos'] = number_column
#     return df

# def add_aa_position_to_df(df: pd.DataFrame) -> pd.DataFrame:
#     """This function adds the amino acid position to the dataframe, under new column 'aa_pos'
#     Method:
#     Add a new column 'aa_pos' to the dataframe
#     Iterate through the "pos_hg19" column, and for each increasment of 3, add the value to the new column
#     """
#     # Initialize variables to keep track of the last value and the incremental number
#     last_value = None
#     incremental_number = 0
#     number_column = []
#
#     # Iterate through the pos_hg19 column and calculate the numbers column
#     for value in df['pos_hg19']:
#         if last_value is None:
#             last_value = value
#             incremental_number = 1
#         elif value - last_value >= 3:
#             incremental_number += 1
#             last_value = value
#         number_column.append(incremental_number)
#
#     # Add the numbers column to the DataFrame
#     df['aa_pos'] = number_column
#     return df


def add_aa_position_to_df_according_to_seq(df: pd.DataFrame, seq: str) -> pd.DataFrame:
    subseq = seq
    positions_list = list(range(1, len(seq)+1)) #This list is consists of positions corresponding with the aa positions in the sequence.
    output_df = pd.DataFrame()
    subdf = df

    print(f"seq is {seq}")
    for aa in seq:
        pos = positions_list[seq.index(aa)]
        print(f"pos is {pos}, aa is {aa}")
        if not is_next_aa_the_same(subseq):
            intermidiate_df = add_position_to_df_by_aa(subdf, pos, aa)
            # add intermidiate_df to output_df
            output_df = pd.concat([output_df, intermidiate_df], ignore_index=True)
            subdf = output_df[output_df['position'].isnull()]
            subseq = subseq[1:]
        if is_next_aa_the_same(subseq):
            # I will only add by the DNA, and will not treat the next aa that is the same! It will loop again for this.
            intermidiate_df = add_position_to_df_when_multiple_concequitive_aa(subdf, pos) # Intermidiate_df will have only the rows with position.
            output_df = pd.concat([output_df, intermidiate_df], ignore_index=True)
            subdf = output_df[output_df['position'].isnull()]
            subseq = subseq[1:]
    df = output_df
    return df


def add_position_to_df_when_multiple_concequitive_aa(subdf: pd.DataFrame, pos) -> pd.DataFrame:
    """
    This method adds the pos to the aa by the addition in the dna column.
    This will only add pos in position column for the aa in this pos.
    """
    last_dna_value = None
    for index, row in subdf.iterrows():
        value = row['pos_hg19']
        if last_dna_value is None:
            last_dna_value = row['pos_hg19']
        elif value - last_dna_value >= 3:
            last_dna_value = value
            row["position"] = pos
            break
        row["position"] = pos
    return subdf


def add_position_to_df_by_aa(subdf: pd.DataFrame, pos, aa):
    for index, row in subdf.iterrows():
        if row["aa_wt"] == aa:
            row["position"] = pos
        else:
            break
    return subdf


def is_next_aa_the_same(seq: str) -> bool:
    if len(seq) < 2:
        return False
    elif seq[0] == seq[1]:
        return True
    return False


def get_consequitive_num(subseq: str) -> int:
    aa = subseq[0]
    count = 0
    for i in subseq:
        if subseq[i] != aa:
            break
        else:
            count+=1
    return count

#
# def test_fit_to_seq(df, seq):
#     for pos, letter in enumerate(seq, 1):
#         # Go over the position and aa_wt in every row and compare to pos and letter.
#         matching_rows = df.loc[df["aa_pos"] == pos]
#
#         # Check if the 'aa_wt' column in the matching rows contains 'letter'
#         if not matching_rows["aa_wt"].isin([letter]).all():
#             print(f"ERROR - pos {matching_rows['aa_pos'].iloc[0]} letter {letter}\n")
#         print(f"{pos} {letter}")


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
    revel_data.head()
    revel_data.columns = ['chr', 'pos_hg19', 'pos_grch38', 'ref_na', 'alt_na', 'aa_wt', 'aa_mut', 'revel_score', "transcript_id"]
    revel_data.head(20)
    seq = uni.get_sequence(gene)

    if reverse:
        revel_data = revel_data[::-1].reset_index(drop=True)
        revel_data = add_aa_position_to_df_reverse(revel_data)
    else:
        revel_data = add_aa_position_to_df_according_to_seq(revel_data, seq)

    # test_fit_to_seq(revel_data, seq)
    # test_fit_to_seq(seq)
    # revel_data = create_variant_column(revel_data)
    revel_data.to_csv(f'{gene}_revel_with_pos.csv', index=False)