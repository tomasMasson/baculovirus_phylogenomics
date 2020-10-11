#!/usr/bin/env python3

import argparse
import subprocess
import shutil
from pathlib import Path
from Bio import SeqIO


def extract_sequences(groups_file, multifasta):
    """
    Extract the fasta sequences for each group and
    write them into a single file.
    """

    # Create a temporal directory to store sequences
    directory = Path.cwd() / "temporal"
    # Store sequence groups as a dictionary
    with open(groups_file, "r") as fh:
        groups = {line.split(",")[0]: line.strip().split(",")[1:]
                  for line in fh
                  if len(line.split(",")) > 2}

    # Store fasta sequences as a dictionary
    seqs = SeqIO.to_dict(SeqIO.parse(multifasta, 'fasta'))
    # For each protein group
    for group in groups:
        # Create a separate file
        filename = directory / group
        with open(filename, "w") as fh:
            # Add sequences based on identifiers
            for identifier in groups[group]:
                seq = str(seqs[identifier].seq)
                fh.write(f">{identifier}\n{seq}\n")


def build_hmm_model(seqs):
    """
    Builds a HMM models for a multifasta file.
    """

    mafft_align = f'mafft --localpair --maxiterate 1000 {seqs} >> {seqs}.aln'
    build_hmm = f'hmmbuild {seqs}.hmm {seqs}.aln'
    subprocess.run(mafft_align, shell=True)
    subprocess.run(build_hmm, shell=True)


def build_hmm_database(db_name):
    '''
    Takes a CSV with sequence groups, extract the sequences from
    a multifasta, performs a multiple sequence alignment and
    prints a HMM database containing all models derived from the alignments.
    '''

    # Iterate over files in the temporal directory
    directory = Path.cwd() / "temporal"
    files = Path.iterdir(directory)
    for file in files:
        # Build orthogroup alignments and HMM models
        build_hmm_model(str(file))
    # Merge all models into a single HMM database and compress it
    hmm_models = str(directory / "*.hmm")
    concatenate = f'cat {hmm_models} > {db_name}'
    compress = f'hmmpress {db_name}'
#    export = f'mv hmm_db results/'
    subprocess.run(concatenate, shell=True)
    subprocess.run(compress, shell=True)
#    subprocess.run(export, shell=True)


def main():
    parser = argparse.ArgumentParser(
        description='Takes a CSV file with sequences groups and creates a HMM database containing all models derived from the multiple sequence alignments',
        usage='build_hmm_database <groups> <multifasta> -o <database>')
    parser.add_argument('groups', help='CSV sequence groups')
    parser.add_argument('multifasta', help='Fasta sequences file')
    parser.add_argument('-o',
                        '--output',
                        nargs='?',
                        default='hmm_db',
                        help='Output file name')
    args = parser.parse_args()

    # Check if a temporal directory exist
    directory = Path.cwd() / "temporal"
    if Path.exists(directory):
        # Delete previous temporal directory
        shutil.rmtree(directory)
    else:
        # Create temporal directory
        Path.mkdir(directory)

    # Extract fasta sequences groups
    extract_sequences(args.groups, args.multifasta)
    # Build the HMM database
    build_hmm_database(args.output)
    # Clean temporal directory
    shutil.rmtree(directory)


if __name__ == '__main__':
    main()
