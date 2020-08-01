#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import os
import subprocess


def sequences_number(seqfile, n_seq):
    """
    Return the number of sequences in a fasta file.
    """
    sequences = list(SeqIO.parse(seqfile, 'fasta'))
    return len(sequences) >= n_seq


def build_hmm_model(seqs):
    """
    Builds a HMM models for a multifasta file.
    """
    mafft_align = f'mafft --localpair --maxiterate 1000 {seqs} >> {seqs[:-3]}.aln'
    build_hmm = f'hmmbuild {seqs[:-3]}.hmm {seqs[:-3]}.aln'
    subprocess.run(mafft_align, shell=True)
    subprocess.run(build_hmm, shell=True)


def create_hmm_database(files_dir, min_number, db_name):
    '''
    Takes a directory path containing several multifasta files
    and prints a HMM database containing all models derived
    from sequence alignments.
    The minimal length is the minimum number of sequences required to
    process the multifasta.
    '''
    work_dir = os.path.abspath(files_dir)
    files = os.listdir(files_dir)
    # Build orthogroup alignments and HMM models
    for file in files:
        seqfile = os.path.join(work_dir, file)
        if seqfile.endswith('.fa') and sequences_number(seqfile, min_number):
            build_hmm_model(seqfile)
#            command1 = f'mafft --localpair --maxiterate 1000 {seqfile} >> {seqfile[:-3]}.aln'
#            command2 = f'hmmbuild {seqfile[:-3]}.hmm {seqfile[:-3]}.aln'
#            subprocess.run(command1, shell=True)
#            subprocess.run(command2, shell=True)
    # Merge all models into a single HMM database and compress it
    concatenate = f'cat {work_dir}/*.hmm > {db_name}'
    compress = f'hmmpress {db_name}'
    export = f'mv {db_name}* results/'
    subprocess.run(concatenate, shell=True)
    subprocess.run(compress, shell=True)
    subprocess.run(export, shell=True)


def main():
    parser = argparse.ArgumentParser(
    description='Takes a directory path containing several multifasta files and creates a HMM database containing all models derived from sequence alignments',
    usage='create_hmm_database <files_dir> <min_number> <output>')
    parser.add_argument('files_dir', help='Directory with multifasta files')
    parser.add_argument('number', type=int, help='Minimal number of sequences per file')
    parser.add_argument('output', help='Output file name')
    args = parser.parse_args()
    create_hmm_database(args.files_dir, args.number, args.output)


if __name__ == '__main__':
    main()
