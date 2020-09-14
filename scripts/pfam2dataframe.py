#!/usr/bin/env python3

import argparse
from Bio import SearchIO


def pfam2dataframe(filename):
    '''
    Returns the first hit of each query in a HMM search.
    Output columns format:
    Protein,Domain,Pfam_ID,E-value,Domain_E-value,Sequence_coverage
    '''
    search = SearchIO.parse(filename, 'hmmscan3-domtab')
    print("Protein,Domain,Pfam_ID,E-value,Domain_E-value,Sequence_coverage")
    for query in search:
        identifier = "_".join(query.id.split("_")[3:5])
        try:
            domain_name = query[0].id
            domain_id = query[0].accession
            evalue = query[0].evalue
            hsp = query[0][0]
            evalue_dom = hsp.evalue
            coverage = (hsp.query_end - hsp.query_start) / query.seq_len * 100
            print(f'{identifier},{domain_name},{domain_id},{evalue},{evalue_dom},{coverage}')
        except IndexError:
            print(f'{identifier},NaN,NaN,NaN,NaN')


def main():
    parser = argparse.ArgumentParser(
            description='Return the first hit for each query of a HMM search',
            usage='pfam2dataframe.py <hmm_search_results>',
            epilog="""
            """)
    parser.add_argument('hmmfile')
    args = parser.parse_args()
    pfam2dataframe(args.hmmfile)


if __name__ == '__main__':
    main()
