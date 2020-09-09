#!/usr/bin/env python3

from Bio import SearchIO
import argparse


def extract_gene_clusters(hmmer_search, reference_names):
    '''
    Complete orthogroups using the results from a hhmscan search
    against a HMM profiles database created with OrthoFinder
    orthogroups.
    '''

    # Load hhmsearch results
    search = list(SearchIO.parse(hmmer_search, 'hmmer3-tab'))
    # Initialize a dict to store the results
    orthogroups = {}
    # Iterate over all query proteins
    for query in search:
        # Retrieve the best hit (Orthogroup) in the database for the query
        try:
            hit_id = query[0].id
            # If the orthogroup doesn't exist , create it and store results
            if hit_id not in orthogroups.keys():
                orthogroups[hit_id] = {'protein_id': [query.id]}
            # If the orthogroup exists, only append results
            elif hit_id in orthogroups.keys():
                orthogroups[hit_id]['protein_id'].append(query.id)
        # Handle queries without hits
        except IndexError:
            pass

    # Create a dict containing reference gene names
    with open(reference_names, "r") as f:
        next(f)
        gene_names = {line.split(",")[0]: line.rstrip().split(",")[1]
                      for line in f}

    # Rename orthogroups using the reference gene names
    for key in orthogroups:
        for protein in orthogroups[key]["protein_id"]:
            # Check if any protein inside the orthogroup has a reference name
            protein_id = "_".join(protein.split("_")[3:5])
            if protein_id in gene_names.keys():
                orthogroups[key]["gene"] = gene_names[protein_id]
                break
            else:
                orthogroups[key]["gene"] = None

    # Results to standard output
    print("Protein,Specie,Orthogroup,Gene_Cluster")
    for key in orthogroups:
        for protein in orthogroups[key]["protein_id"]:
            # Extract protein and specie identifier
            protein_id = "_".join(protein.split("_")[3:5])
            specie = "NC_" + "_".join(protein.split("_")[1:2])
            gene_cluster = orthogroups[key]["gene"]
            print(f"{protein_id},{specie},{key},{gene_cluster}")


def main():
    parser = argparse.ArgumentParser(
    description='Returns to standard output the proteins contained in each orthogroup.',
    usage='parse_hmmer-tab <hmmer> <names_file>')
    parser.add_argument('hmmer', type=str, help='HMMER search in tab format')
    parser.add_argument('names', type=str, help='File containig sequence names')
    args = parser.parse_args()
    extract_gene_clusters(args.hmmer, args.names)


if __name__ == '__main__':
    main()
