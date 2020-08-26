#!/usr/bin/env python3

from Bio import SearchIO
import argparse


def add_missing_orthologs(hmmer_search):
    '''
    Complete orthogroups using the results from a hhscan search
    against a HMM profiles database created with OrthoFinder
    orthogroups.
    '''

    # Load hhmsearch results
    search = list(SearchIO.parse(hmmer_search, 'hmmer3-tab'))
    # Initialize a dict to store the results
    orthogroups = {}
    # Iterate over all query proteins
    for query in search:
        # Store query specie and protein name
        query_protein = '_'.join(query.id.split('_')[3:5])
        query_specie = 'NC_' + str(query.id.split('_')[1])
        # Retrieve the best hit (Orthogroup) in the database for the query
        try:
            hit_id = query[0].id
            # If the orthogroup doesn't exist , create it and store results
            if hit_id not in orthogroups.keys():
                orthogroups[hit_id] = {
                        'species': [query_specie],
                        'proteins': [query_protein],
                        'full_ids': [query.id]
                        }
            # If the orthogroup exists, only append results
            elif hit_id in orthogroups.keys():
                orthogroups[hit_id]['species'].append(query_specie)
                orthogroups[hit_id]['proteins'].append(query_protein)
                orthogroups[hit_id]['full_ids'].append(query.id)
        # Handle queries without hits
        except IndexError:
            pass
    # Results dict to standard output
    for key in orthogroups:
        print(f'{key},{",".join(orthogroups[key]["full_ids"])}')
    return orthogroups


def main():
    parser = argparse.ArgumentParser(
    description='Returns to standard output the proteins contained in each orthogroup.',
    usage='add_missing_orthologs <hmmer>')
    parser.add_argument('hmmer', type=str, help='HMMER search in tab format')
    args = parser.parse_args()
    add_missing_orthologs(args.hmmer)


if __name__ == '__main__':
    main()
