#!/usr/bin/env python3

from Bio import SeqIO
import argparse


def get_taxonomic_distribution(multifasta):
    '''
    Returns the a list of unique species present inside an
    orthogroup. Additionally, print the number of sequences per
    specie.
    '''
    sequences = SeqIO.parse(multifasta, 'fasta')
    dic = {}
    for seq in sequences:
        specie_id = 'NC_' + str(seq.id.split('_')[1])
        protein_id = '_'.join(seq.id.split('_')[3:5])
        dic[specie_id] = protein_id

    if 'NC_001623.1' in dic.keys():
        orthogroup = dic['NC_001623.1']
    elif 'NC_002816.1' in dic.keys():
        orthogroup = dic['NC_002816.1']
    else:
        orthogroup = multifasta.split('.fa')[0]
    n_total_proteins = len(dic)
    n_unique_proteins = len(set(dic.keys()))
    n_unique_species = len(set(dic.values()))
    print(orthogroup, n_total_proteins, n_unique_proteins, n_unique_species)


def check_missing_species(multifasta, reference_file):
    '''
    Returns the species present in the reference file but not
    represented inside the multifasta.
    '''
    with open(reference_file, 'r') as f:
        references = [line.strip() for line in f]

    sequences = SeqIO.parse(multifasta, 'fasta')
    species = ['NC_'+str(seq.id.split('_')[1]) for seq in sequences]

    print(set(references).difference(set(species)))


def main():
    parser = argparse.ArgumentParser(
    description='Returns the number of unique species represented inside a multifasta and the missing species qhen compared against a reference file',
    usage='taxonomy_distribution <multifasta> <reference>')
    parser.add_argument('multifasta', type=str, help='Multifasta file')
    parser.add_argument('references', type=str, help='Reference species file, with one record per line')
    args = parser.parse_args()
    get_taxonomic_distribution(args.multifasta)
    check_missing_species(args.multifasta, args.references)


if __name__ == '__main__':
    main()
