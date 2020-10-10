#!/usr/bin/env python3

import argparse
import sys
import pandas as pd


def extract_pfam_groups(interpro_search):
    """
    Extract the Pfam instances present in a InterPro
    results file (TSV) and groups the sequences that
    shares the same Pfam domain.
    Output goes to the file 'pfam_groups.csv'
    """

    # Define column names
    column_names = [
            "Protein", "Hash", "Length", "DB",
            "DB_instance", "Name", "Start", "End",
            "E-value", "Date", "Description"
            ]
    # Create DataFrame with InterPro results
    df = pd.read_csv(interpro_search,
                     sep="\t",
                     names=column_names)
    # Keep only Pfam instances
    pfam = df[df.loc[:, "DB"] == "Pfam"]
    # Group using Pfam identifier
    groups = pfam.groupby("DB_instance")
    # Recover protein groups
    protein_groups = {label: group.Protein.to_list()
                      for label, group in groups}
    for label in protein_groups:
        proteins = ",".join(protein_groups[label])
        print(f"{label},{proteins}")


def main():
    """
    Command line argument parser
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("interpro_results",
                        help="InterPro search results")
#    parser.add_argument("-o",
#                        "--outfile",
#                        nargs='?',
#                        type=argparse.FileType('w'),
#                        default=sys.stdout)
    args = parser.parse_args()
    extract_pfam_groups(args.interpro_results)


if __name__ == "__main__":
    main()
