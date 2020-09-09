#!/usr/bin/env python3

import argparse
from datetime import datetime
import os
import shutil
import snakemake

# Pipeline general configurations
#date = datetime.now().strftime('%h') + datetime.now().strftime('%d')
snakefile = 'Snakefile'
#orthofinder_dir = os.path.join(f'OrthoFinder/Results_{date}', 'Orthogroup_Sequences')
#results_dir = 'results/'


def main():
#    parser = argparse.ArgumentParser(
#            description='Protein orthology inference pipeline based on Orthofinder and Broccoli',
#            usage='command.py <proteomes> <min_number>',
#            epilog=
#            '''
#            A csv file with gene cluster is created at results folder.
#            Each row displays the best cluster identifier and the individual
#            protein identifiers.
#            ''')
##    parser.add_argument('proteomes', help='Directory containing proteome files')
#    parser.add_argument('min_number', help='Minimal number of sequences required to process a multifasta orthogroup')
#    parser.add_argument('genomes', help='Folder with genome sequences')
#    parser.add_argument('proteomes', help='Folder with proteomes sequences')
#    args = parser.parse_args()
#    genomes = args.genomes
#    proteomes = args.proteomes
#    min_number = args.min_number
#
#    config = {
#            'genomes': genomes,
#            'proteomes': proteomes,
#            'orthof_dir': orthofinder_dir,
#            'min_number': min_number
#            }
    if os.path.isdir('results'):
        shutil.rmtree('results')
        print('Results directory cleaned')
    else:
        os.mkdir('results')
    analysis = snakemake.snakemake(snakefile,
                                   printshellcmds=True
                                   )
    if analysis:
        return "Analysis Completed"
    return "Something went wrong"


if __name__ == '__main__':
    main()
