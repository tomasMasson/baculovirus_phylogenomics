#!/usr/bin/env python3

from Bio import SeqIO
import argparse

reference_species = {}

core_genes = {
        'NP_054035.1': 'lef2',
        'NP_054043.1': 'lef1',
        'NP_054051.1': 'pif1',
        'NP_054069.1': 'p47',
        'NP_054079.1': 'lef8',
        'NP_054082.1': 'ac53',
        'NP_054084.1': 'vp1054',
        'NP_054092.1': 'lef9',
        'NP_054095.1': 'dna-pol',
        'NP_054096.1': 'desmop',
        'NP_054098.1': 'ac68',
        'NP_054107.1': 'vlf1',
        'NP_054108.1': 'ac78',
        'NP_054110.1': 'gp41',
        'NP_054111.1': 'ac81',
        'NP_054113.1': 'vp91',
        'NP_054119.1': 'vp39',
        'NP_054120.1': 'lef4',
        'NP_054122.1': 'p33',
        'NP_054123.1': 'p18',
        'NP_054124.1': 'odv-e25',
        'NP_054125.1': 'helicase',
        'NP_054126.1': 'ac96',
        'NP_054128.1': '38k',
        'NP_054129.1': 'lef5',
        'NP_054130.1': 'p6.9',
        'NP_054131.1': 'p40',
        'NP_054133.1': 'p48',
        'NP_054139.1': 'odv-ec43',
        'NP_054140.1': 'ac110',
        'NP_054145.1': 'pif3',
        'NP_054149.1': 'pif1',
        'NP_054163.1': 'alk-exo',
        'NP_054168.1': 'p74',
        'NP_054173.1': '49k',
        'NP_054174.1': 'odv-e18',
        'NP_054175.1': 'odv-e27',
        'NP_054179.1': 'pif5',
        }


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
#    if 'NC_001623.1' in dic.keys():
#        orthogroup = dic['NC_001623.1']
#    elif 'NC_002816.1' in dic.keys():
#        orthogroup = dic['NC_002816.1']
    for key in core_genes:
        if key in dic.values():
            orthogroup = core_genes[key]
            break
        else:
            orthogroup = 'NA'
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

    print(','.join(set(references).difference(set(species))))


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
