#!/usr/bin/env python3

import argparse
from datetime import datetime
import os
import shutil
import snakemake

date = datetime.now().strftime('%h') + datetime.now().strftime('%d')
snakefile = 'Snakefile'
og_dir = os.path.join(f'Results_{date}', 'Orthogroup_Sequences')


def main():
    parser = argparse.ArgumentParser(
            description='Protein orthology inference pipeline based on Orthofinder and Broccoli',
            usage='command.py <proteomes>')
    parser.add_argument('proteomes',
                        help='Directory containing proteomes fasta files. Extension of the file must be .fasta (Broccoli requirement)')
    args = parser.parse_args()
    data_dir = args.proteomes

    config = {
            'data_dir': data_dir,
            'og_dir': og_dir
            }
    if os.path.isdir('results'):
        shutil.rmtree('results')
        print('Results directory cleaned')
    else:
        os.mkdir('results')
    analysis = snakemake.snakemake(snakefile, printshellcmds=True,
                                   config=config
                                   )
    if analysis:
        return 0
    return 1


if __name__ == '__main__':
    main()
