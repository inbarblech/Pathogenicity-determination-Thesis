from data_retrievel_and_feature_extraction import uniprot_info as uni
import pandas as pd

def transofrm_go_terms_to_df(go_terms: dict) -> pd.DataFrame:
    data = pd.DataFrame.from_dict(go_terms, orient='index')
    data = data.reset_index()
    data = data.drop(columns=['index'])
    return data


def split_go_terms_to_categories(df: pd.DataFrame) -> tuple:
    """Split the GO terms into 3 dataframes, according to the 3 GO categories, and save them to csv.
    Args:
        df (pd.DataFrame): A dataframe containing the GO terms for the given UniProt accession ID
    Returns:
        tuple: A tuple containing the 3 dataframes, each containing the GO terms for the specific gene in the df.
        first dataframe: cellular component
        second dataframe: biological process
        third dataframe: molecular function
    """

    df_cellular_component = df[df[0].str.startswith('C:')]
    df_biological_process = df[df[0].str.startswith('P:')]
    df_molecular_function = df[df[0].str.startswith('F:')]

    return df_cellular_component, df_biological_process, df_molecular_function


if __name__ == "__main__":
    go_terms = uni.get_go_terms('BRCA1')
    go_terms_df = transofrm_go_terms_to_df(go_terms)
    df_c, df_p, df_f = split_go_terms_to_categories(go_terms_df)
    # save go terms to csv
    df_c.to_csv('cellular_component.csv')
