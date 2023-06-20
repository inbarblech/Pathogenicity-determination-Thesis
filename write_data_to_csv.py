import pandas as pd


def create_new_df():
    df = pd.DataFrame(columns=['gene', 'pathogenicity', 'uniprot_id', 'variant', 'pathogenicity', 'delta_stability',
                               'delta_energy'])
    return df


def add_variant_data_to_df(df, variant_data: dict):
    """Adds the variant data to the csv file."""
    df = df.append(variant_data, ignore_index=True)
    return df
