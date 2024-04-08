import pandas as pd
import requests
import json
import sys


def get_genomic_location(gene_symbol, species="homo_sapiens"):
    url = f"https://rest.ensembl.org/lookup/symbol/{species}/{gene_symbol}?expand=1"

    # Send a GET request to the Ensembl API
    response = requests.get(url, headers={"Content-Type": "application/json"})

    if response.ok:
        # Parse the JSON response
        data = response.json()

        # Extract chromosome, start, end, and strand information
        chromosome = data['seq_region_name']
        start = data['start']
        end = data['end']
        strand = data['strand']
        assembly_name = data['assembly_name']

        return chromosome, start, end, strand, assembly_name
    else:
        print(f"Failed to retrieve genomic location for {gene_symbol}. Status code: {response.status_code}")
        return None


if __name__ == "__main__":
    genes = ["COL2A1", "COL4A5"]

    for gene in genes:
        # Get locations
        chromosome, start, end, strand_num, assmebly_name = get_genomic_location(gene)
        print(start)

