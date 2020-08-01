#!/usr/bin/env python3

import entrezpy.conduit
import argparse
import sys


def query_entrez(accession, rettype):
    '''
    Retrieves the full sequence (genome/proteome) for the accession number
    provided.
    Output goes to standard output.
    '''
    # Replace email with yours
    email = 'tomasmasson0@gmail.com'
    # Create a Conduit instance
    search = entrezpy.conduit.Conduit(email)
    # Create a new search pipeline
    fetch_data = search.new_pipeline()
    # Set Esearch query
    sid = fetch_data.add_search({'db': 'nucleotide',
                                 'term': accession}
                                )
    # Set Efetch query
    fetch_data.add_fetch({'retmax': 10, 'retmode': 'text',
                          'rettype': rettype},
                          dependency=sid
                         )
    # Execute search and write results to output file
    with open(f'{accession}.fa', 'w') as f:
        sys.stdout = f
        search.run(fetch_data)


def main():
    parser = argparse.ArgumentParser(
            description='Download NCBI sequences based on an accession numbers list',
            usage='retrieve_sequence.py <acc_list> <rettype>')
    parser.add_argument('accessions', help='Accession numbers list')
    parser.add_argument('rettype', help='Type of data retrieved. Options available are the same as efetch from NCBI E-Utilities')
    args = parser.parse_args()
    accession_list = args.accessions
    rettype = args.rettype

    with open(accession_list, 'r') as fh:
        accessions = [line.strip() for line in fh]
    for accession in accessions:
        query_entrez(accession, rettype)


if __name__ == '__main__':
    main()
