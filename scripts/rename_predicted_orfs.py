#!/usr/bin/env python3

import argparse
import subprocess
from Bio import SearchIO
from Bio import SeqIO


def diamond_search(query, database):
    """
    Builds a BLAST database and performs a BLASTp search
    againts query sequences.
    """

    build_db = f'diamond makedb --in {database} --db tmp.db'
    run_diamond = f'diamond blastp --threads 4 --query {query} --db tmp.db --outfmt 5 --out tmp.xml'
    subprocess.run(build_db, shell=True)
    subprocess.run(run_diamond, shell=True)


def create_names_dict(blast_search):
    """
    Takes the results from a BLASTp search and return the best hit id for
    each query. If there are no hits, the query id is return.
    """

    search = SearchIO.parse(blast_search, 'blast-xml')
    names_dict = {}
    for query in search:
        if query.hits:
            names_dict[query.id] = query.hits[0].id
        else:
            names_dict[query.id] = query.id
    return names_dict


def rename_orfs(query, database, output):
    """
    """
    diamond_search(query, database)
    names_dict = create_names_dict("tmp.xml")
    seqs = SeqIO.parse(query, "fasta")
    with open(output, "w") as f:
        for seq in seqs:
            if seq.id in names_dict.keys():
                f.write(f">{names_dict[seq.id]}\n{seq.seq}\n")


def main():
    parser = argparse.ArgumentParser(
    description='Performs a local BLASTp search using the query and database given',
    usage='renama_predicted_orfs <query> <database> <output>')

    parser.add_argument('query', help='Query sequences in fasta format')
    parser.add_argument('database', help='Sequence to build the database')
    parser.add_argument('output', help='Output file name')
    args = parser.parse_args()
    rename_orfs(args.query, args.database, args.output)
    subprocess.run("rm -rf tmp.db* tmp.xml", shell=True)


if __name__ == '__main__':
    main()
