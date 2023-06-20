import requests as req
import re


def get_uniprot_url(gene_name):
    """Returns the URL for the Uniprot page for the given gene name."""
    url = f"https://rest.uniprot.org/uniprotkb/search?query=(gene:{gene_name})%20AND%20(taxonomy_id:9606)%20AND%20(reviewed:true)"
    return url


def get_uniprot_id(gene_name):
    """Returns the Uniprot ID for the given gene name."""
    data = get_uniprot_json(gene_name)
    primary_accession = data['results'][0]['primaryAccession']
    return primary_accession


def get_uniprot_json(gene_name):
    # The base URL for UniProt's search API
    url = get_uniprot_url(gene_name)
    # Make a request to the search API
    response = req.get(url)
    # Extract the JSON data from the response
    data = response.json()
    return data


def get_sequence(gene_name):
    """Returns the sequence for the given Uniprot"""
    data = get_uniprot_json(gene_name)
    sequence = data['results'][0]['sequence']['value']
    return sequence


def get_sequence_length(gene_name):
    """Returns the length of the sequence for the given Uniprot"""
    data = get_uniprot_json(gene_name)
    length = data['results'][0]['sequence']['length']
    return length


