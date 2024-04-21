import json

import pandas as pd
import requests as req
import re
import csv
from io import StringIO


def get_uniprot_url(gene_name) -> str:
    """Returns the URL for the Uniprot page for the given gene name."""
    url = f"https://rest.uniprot.org/uniprotkb/search?query=(gene:{gene_name})%20AND%20(taxonomy_id:9606)%20AND%20(reviewed:true)"
    return url


def get_gene_name_from_uniprot_id(uniprot_id) -> str:
    """Returns the gene name for the given Uniprot ID."""
    data = get_uniprot_json(uniprot_id)
    gene_name = data['results'][0]['gene']
    return gene_name


def get_uniprot_id(gene_name) -> str:
    """Returns the Uniprot ID for the given gene name."""
    data = get_uniprot_json(gene_name)
    primary_accession = data['results'][0]['primaryAccession']
    return primary_accession


def get_uniprot_json(gene_name) -> dict:
    # The base URL for UniProt's search API
    url = get_uniprot_url(gene_name)
    # Make a request to the search API
    response = req.get(url)
    # Extract the JSON data from the response
    data = response.json()
    return data


def write_uniprot_json_to_file(gene_name):
    # The base URL for UniProt's search API
    url = get_uniprot_url(gene_name)
    # Make a request to the search API
    response = req.get(url)
    # Extract the JSON data from the response
    data = response.json()

    with open('C:\\Users\\InbarBlech\\Downloads\\data_json.txt', 'w') as f:
        json.dump(data, f)


def get_sequence(gene_name) -> str:
    """Returns the sequence for the given Uniprot"""
    data = get_uniprot_json(gene_name)
    sequence = data['results'][0]['sequence']['value']
    return sequence


def get_type_of_residue_membrane_or_globular(gene_name, position) -> bool:
    """Returns True if residue is in transmembranal area, False if in globular."""
    data = get_uniprot_json(gene_name)
    cellular_location_info = data['results'][0]['features']
    for feature in cellular_location_info:
        if feature['type'] == "Transmembrane":
            if feature['location']['start']['value'] <= position <= feature['location']['end']['value']:
                return True
    return False


def is_protein_membranal(gene_name) -> bool:
    """Returns True if there is any transmembrane part to the protein.
    False otherwise."""
    data = get_uniprot_json(gene_name)
    cellular_location_info = data['results'][0]['features']
    for feature in cellular_location_info:
        if feature['type'] == "Transmembrane":
            return True
    return False


# def get_secondary_structure(gene_name, position):
#     """Returns the secondary structure of the protein in the given position."""
#     data = get_uniprot_json(gene_name)
#     secondary_structure = data['results'][0]['features']
#     for feature in secondary_structure:
#         if feature['type'] == 'Helix' or feature['type'] == 'Beta strand' or feature['type'] == 'Turn':
#             if feature['location']['start']['value'] <= position <= feature['location']['end']['value']:
#                 return feature['type']
#     return 'Loop'


def get_go_terms(gene_name) -> dict:
    """Returns the GO terms for the given UniProt accession ID
    in a dictionary where the keys are the GO IDs and the values are the GO terms."""
    data = get_uniprot_json(gene_name)
    go_terms = dict()

    for entry in data['results'][0]["uniProtKBCrossReferences"]:
        if entry["database"] == "GO":
            go_id = entry["id"]
            for property in entry["properties"]:
                if property["key"] == "GoTerm":
                    go_term = property["value"]
                    if go_id in go_terms:
                        go_terms[go_id].append(go_term)
                    else:
                        go_terms[go_id] = [go_term]
    return go_terms


# def export_go_terms_to_csv(data, output_file):
#     with open(output_file, 'w', newline='') as file:
#         writer = csv.writer(file)
#         writer.writerow(['ID', 'Category', 'Property'])
#
#         for go_id, properties in data.items():
#             for property in properties:
#                 category = property[0]
#                 property_value = property[2:]
#                 writer.writerow([go_id, category, property_value])
#
#     print(f"Data has been exported to '{output_file}'.")


def get_go_terms_csv(data : dict) -> str:
    """Returns a CSV string containing the GO terms for the given UniProt accession ID
    in the format: ID, Category, Property
    Args:
        data (dict): A dictionary containing the GO terms for the given UniProt accession ID
    Returns:
        str: A CSV string containing the GO terms for the given UniProt accession ID

    Call this function in the following way:
    csv_data = get_go_terms_csv(go_terms)
    """
    output = StringIO()
    writer = csv.writer(output)
    writer.writerow(['ID', 'Category', 'Property'])

    for go_id, properties in data.items():
        for property in properties:
            category = property[0]  # Extract the first character as the category
            property_value = property[2:]  # Exclude the first two characters (e.g., "C:", "P:", "F:")
            writer.writerow([go_id, category, property_value])

    output.seek(0)  # Move the stream cursor to the beginning of the output
    return output.getvalue()


def get_sequence_length(gene_name) -> int:
    """Returns the length of the protein for the given Uniprot, using the sequence length from Uniprot."""
    data = get_uniprot_json(gene_name)
    length = data['results'][0]['sequence']['length']
    return length


def get_secondary_structure(gene_name, position):
    """Returns the secondary structure of the protein in the given position."""
    data = get_uniprot_json(gene_name)
    secondary_structure = data['results'][0]['features']
    for feature in secondary_structure:
        if feature['type'] == 'Helix' or feature['type'] == 'Beta strand' or feature['type'] == 'Turn':
            if feature['location']['start']['value'] <= position <= feature['location']['end']['value']:
                return feature['type']
    return 'Loop'


if __name__ == "__main__":
    go_terms = get_go_terms("BRCA1")

    # save dict as df
    df = pd.DataFrame.from_dict(go_terms, orient='index')

    # split df into 3 dataframes, according to the 3 go categories, and save them to csv
    df_cellular_component = df[df[0].str.startswith('C:')]
    df_biological_process = df[df[0].str.startswith('P:')]
    df_molecular_function = df[df[0].str.startswith('F:')]



    # # save go terms to csv
    # csv_data = get_go_terms_csv(go_terms)
    # # save csv to file
    # with open("go_terms.csv", "w") as file:
    #     file.write(csv_data)
