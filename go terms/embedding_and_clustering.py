import pandas as pd
from sentence_transformers import SentenceTransformer
import numpy as np


# Preprocessing functions
def delete_sign_from_column(dataframe, column_name, sign_to_remove):
    dataframe[column_name] = dataframe[column_name].str.replace(sign_to_remove, "")
    return dataframe


def split_to_list(dataframe, column_name):
    dataframe[column_name] = dataframe[column_name].str.split(', ')
    return dataframe


# Embedding functions
def get_embeddings(term):
    if term is ['']:
        return None
    # Get the embeddings for the term
    embedding = model.encode(term)
    return embedding


def get_embeddings_for_cell(cell):
    if cell is ['']:
        return None
    mean_embeddings = np.mean([get_embeddings(term) for term in cell], axis=0)
    return mean_embeddings


def write_embeddings_for_all_columns(dataframe, column_name):
    dataframe[f"{column_name}_Embedding"] = dataframe[column_name].apply(lambda x: get_embeddings_for_cell(x))
    return dataframe


def get_averaged_embeddings(terms):
    if not terms:
        return None
    # Join the terms into a single string
    text = ' '.join(terms)
    # Get the embeddings for the combined string
    embedding = model.encode(text)
    averaged_embedding = np.mean(embedding, axis=0)
    return averaged_embedding


if __name__ == "__main__":
    # Load a pre-trained sentence transformer model
    model_name = 'bert-base-nli-mean-tokens'  # You can choose a different model
    bio_model = "dmis-lab/biobert-v1.1"
    model = SentenceTransformer(model_name)

    ### Loading and preprocessing ###
    # Load your CSV file
    input_csv_file_path = 'C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\hl_genes_with_go_terms_lists.csv'
    df = pd.read_csv(input_csv_file_path)

    # remove all signs from the go terms, and split to lists
    for column in ['Cellular Component', 'Biological Process', 'Molecular Function']:
        for sign in ["'", "[", "]"]:
            df = delete_sign_from_column(df, column, sign)
        df = split_to_list(df, column)

    ### Embedding ###
    # Apply the embeddings to each row in the DataFrame
    for column in ['Cellular Component', 'Biological Process', 'Molecular Function']:
        df[column] = df[column].apply(lambda x: write_embeddings_for_all_columns(df, column))

    # Save the DataFrame to a CSV file
    csv_file_path = f'C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\go terms\\hl_genes_with_go_terms_mean_embeddings_29.1.csv'
    df.to_csv(csv_file_path, index=False)