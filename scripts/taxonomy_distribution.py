#!/usr/bin/env python3

from Bio import SeqIO
import argparse


def get_taxonomic_distribution(multifasta):
    '''
    Returns the a list of unique species present inside an
    orthogroup. Additionally, print the number of sequences per
    specie.
    '''
    orthogroup = multifasta.split('.fa')[0]
    sequences = SeqIO.parse(multifasta, 'fasta')
    dic = {}
    for seq in sequences:
        specie_id = 'NC_' + str(seq.id.split('_')[1])
        protein_id = '_'.join(seq.id.split('_')[3:5])
        dic[protein_id] = specie_id
    total_proteins = len(dic)
    unique_proteins = len(set(dic.keys()))
    unique_species = len(set(dic.values()))
    print(orthogroup, total_proteins, unique_proteins, unique_species)

def main():
    parser = argparse.ArgumentParser(
    description='Returns the number of unique species represented inside a multifastaa',
    usage='taxonomy_distribution <multifasta>')
    parser.add_argument('multifasta', type=str, help='Multifasta file')
    args = parser.parse_args()
    get_taxonomic_distribution(args.multifasta)


if __name__ == '__main__':
    main()
