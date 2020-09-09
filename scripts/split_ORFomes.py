#!/usr/bin/env python3

import argparse
from Bio import SeqIO


def split_ORFomes(multifasta):
    """
    Split a multifasta file into individual proteomes
    based on species identifier.
    """
    # Ingest sequences
    seqs = SeqIO.parse(multifasta, "fasta")
    # Initializa proteomes dictionary
    proteomes = {}
    # Store each proteome in the dict
    for seq in seqs:
        specie = "NC_" + seq.id.split("_")[1][:-2]
        if specie not in proteomes:
            proteomes[specie] = [seq]
        elif specie in proteomes:
            proteomes[specie].append(seq)

    # Write each proteome separately to a file
    for proteome in proteomes:
        # Set directory destination
        filename = f"{proteome}.faa"
        with open(filename, "w") as fh:
            for seq in proteomes[proteome]:
                # Fasta format
                fh.write(f">{seq.id}\n{seq.seq}\n")


def main():
    "Command line parser"
    parser = argparse.ArgumentParser(
            description="Splits a multifasta file into individual proteomes",
            usage="split_ORFomes.py <multifasta>")
    parser.add_argument("multifasta", help="File containing all proteins")
    args = parser.parse_args()
    split_ORFomes(args.multifasta)


if __name__ == "__main__":
    main()
