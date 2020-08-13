#!/usr/bin/env python3

import argparse
from Bio import SearchIO


def extract_hmmer_hits(filename):
    '''
    Returns the first hit of each query in a HMM search.
    Output columns format:
    identifier,best_hit,e-value,domain_e-value,coverage
    '''
    search = SearchIO.parse(filename, 'hmmer3-text')
    for query in search:
        identifier = query.id
        try:
            best_hit = query[0].id
            evalue = query[0].evalue
            hsp = query[0][0]
            evalue_dom = hsp.evalue
            coverage = (hsp.query_end - hsp.query_start) / query.seq_len
            print(f'{identifier},{best_hit},{evalue},{evalue_dom},{coverage}')
        except IndexError:
            print(f'{identifier},NA,NA,NA,NA')


def main():
    parser = argparse.ArgumentParser(
            description='Return the first hit for each query of a HMM search',
            usage='extract_hmmer_hits.py <hmm_search_results>',
            epilog="""
            Output columns format:

            identifier,best_hit,e-value,domain_e-value,coverage
            """)
    parser.add_argument('hmmfile')
    args = parser.parse_args()
    extract_hmmer_hits(args.hmmfile)


if __name__ == '__main__':
    main()
