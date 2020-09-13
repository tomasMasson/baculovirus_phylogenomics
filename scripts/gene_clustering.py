#!/usr/bin/env python3

from Bio import SearchIO
import argparse


def gene_clustering(hmmer_search, reference_names):
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

    # Create a dict containing AcMNPV reference gene names
    with open(reference_names, "r") as f:
        next(f)
        reference_gene = {line.split(",")[1]: line.rstrip().split(",")[2]
                          for line in f if line.startswith("acmnpv")}
    # Create a dict containing CpGV gene names
    with open(reference_names, "r") as f:
        next(f)
        cpgv_genes = {line.split(",")[1]: line.rstrip().split(",")[2]
                      for line in f if line.startswith("other")}

    # Rename orthogroups using the reference gene names
    for key in orthogroups:
        for protein in orthogroups[key]["protein_id"]:
            # Check if any protein inside the orthogroup has a reference name
            protein_id = "_".join(protein.split("_")[3:5])
            if protein_id in reference_gene.keys():
                orthogroups[key]["gene"] = reference_gene[protein_id]
                break
            elif protein_id in cpgv_genes.keys():
                orthogroups[key]["gene"] = cpgv_genes[protein_id]
            elif "gene" not in orthogroups[key]:
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
            description='Returns a csv file with the HMMER search results.',
            usage='parse_hmmer-tab <hmmer> <reference_genes>')
    parser.add_argument('hmmer', type=str, help='HMMER search in tab format')
    parser.add_argument('references', type=str, help='Reference names (CSV)')
    args = parser.parse_args()
    gene_clustering(args.hmmer, args.references)


if __name__ == '__main__':
    main()
