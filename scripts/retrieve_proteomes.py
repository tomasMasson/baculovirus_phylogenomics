#!/usr/bin/env python3

import entrezpy.conduit
import argparse
import sys


def retrieve_proteome(taxid):
    '''
    Retrieves the RefSeq proteomes available for
    the taxa specified.
    '''
    # Replace email with yours
    email = 'tomasmasson0@gmail.com'
    # Create a Conduit instance
    search = entrezpy.conduit.Conduit(email)
    # Create a new search pipeline
    fetch_data = search.new_pipeline()
    # Set Esearch query (keep only RefSeq sequences)
    query = f"txid{taxid}[Organism:exp]+refseq[filter]"
    sid = fetch_data.add_search({'db': 'nucleotide',
                                 'term': query})
    # Set Efetch query
    fetch_data.add_fetch({'retmax': 10,
                          'retmode': 'text',
                          'rettype': 'fasta_cds_aa'},
                          dependency=sid)
    # Execute search and return accessions
    peptides = search.run(fetch_data)

    return peptides


def main():
    parser = argparse.ArgumentParser(
        description='Download proteome data for a taxa list',
        usage='retrieve_proteomes.py <taxa_list>')
    parser.add_argument('taxa',
                        help='Taxa list')
    args = parser.parse_args()
    taxa_list = args.taxa
    # Create a taxa list
    with open(taxa_list, 'r') as fh:
        taxa = [taxon.strip() for taxon in fh]
    # Iterate over taxa and fetch proteome data
    with open("proteomes.faa", 'w') as f:
        sys.stdout = f
        for taxon in taxa:
            retrieve_proteome(taxon)


if __name__ == '__main__':
    main()
