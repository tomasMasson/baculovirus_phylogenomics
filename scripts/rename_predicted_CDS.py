#!/usr/bin/env python3

import argparse
from Bio import SeqIO


def annotate_orfs(orfs):
    """
    Takes a multifaste containing ORFs predicted by ORFfinder
    and modify the header to have the same format as NCBI
    proteome fasta.
    0-based index from ORFfinder is changed to 1-based index.

    >lcl|ORF5_NC_XXXXXX.1:25278:25502 -> >lcl|NC_XXXXXX.1_prot_ORF_25279:25503_predicted
    """

    seqs = SeqIO.parse(orfs, "fasta")
    for seq in seqs:
        if "ORF" in seq.id:
            species = "_".join(seq.id.split(":")[0].split("_")[1:])
            start = int(seq.id.split(":")[1]) + 1
            stop = int(seq.id.split(":")[2]) + 1
            new_id = f"lcl|{species}_prot_ORF_{start}:{stop}_predicted"
            print(f">{new_id}\n{seq.seq}")


def main():
    parser = argparse.ArgumentParser(
    description='Takes a fasta file with predicted ORFs from ORFfinder and fixes the headers',
    usage='rename_predicted_orfs <orfs>')
    parser.add_argument('orfs', help='ORFs to rename in fasta format')
    args = parser.parse_args()
    annotate_orfs(args.orfs)


if __name__ == '__main__':
    main()
