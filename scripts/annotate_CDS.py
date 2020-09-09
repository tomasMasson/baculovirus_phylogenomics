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

    # Database build
    build_db = f'diamond makedb --in {database} --db tmp.db'
    # Perform search
    run_diamond = f'diamond blastp --threads 4 --query {query} --db tmp.db --outfmt 5 --out tmp.xml'
    subprocess.run(build_db, shell=True)
    subprocess.run(run_diamond, shell=True)


def retrieve_new_cds(multifasta, blast_search):
    """
    Takes the results from a BLASTp search and return the CDS id for
    the queries without any hit.
    """

    search = SearchIO.parse(blast_search, 'blast-xml')
    # Collect identifiers from BLASTp result
    new_cds_id = []
    for query in search:
        # Only retain sequences without hits
        if not query.hits:
            new_cds_id.append(query.id)
    # Extract sequences from multifasta
    new_cds = {}
    seqs = SeqIO.parse(multifasta, "fasta")
    for seq in seqs:
        # Only sequences not identified in BLASTp
        if seq.id in new_cds_id:
            new_cds[seq.id] = seq.seq
    return new_cds


def update_proteome(cds, proteome, output):
    """
    Update query proteome with predicted ORFs
    with no detectable orthologs in the proteome.
    """
    # Store initial proteome
    seqs = SeqIO.parse(proteome, "fasta")
    # Compute new CDS
    diamond_search(cds, proteome)
    new_cds = retrieve_new_cds(cds, "tmp.xml")
    # Write an updated proteome to output file
    with open(output, "w") as fh:
        for seq in seqs:
            fh.write(f">{seq.id}\n{seq.seq}\n")
        for cds in new_cds:
            fh.write(f">{cds}\n{new_cds[cds]}\n")


def main():
    parser = argparse.ArgumentParser(
    description='Add non-annotated ORFs to a proteome based on BLASTp',
    usage='annotate_predicted_CDS <query> <database> <output>')

    parser.add_argument('query', help='Query sequences in fasta format')
    parser.add_argument('database', help='Sequence to build the database')
    parser.add_argument('output', help='File output name')
    args = parser.parse_args()

    # Run annotation
    update_proteome(args.query, args.database, args.output)
    # Clean intermediate files
    subprocess.run("rm -rf tmp.db* tmp.xml", shell=True)


if __name__ == '__main__':
    main()
