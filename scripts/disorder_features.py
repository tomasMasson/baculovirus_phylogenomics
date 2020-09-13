#!/usr/bin/env python3

import argparse
import pandas as pd


class DisorderFeatures:
    def __init__(self, values):
        """
        Initializes with a list of predicted disorder values
        """
        self.values = values

    def disorder_content(self):
        """
        Filter disordered residues with a predicted
        value greater than 0.5 and compute average
        disorder content
        """
        # Filter residues with score equal or greater than 0.5
        dis_res = [res for res in self.values if float(res) >= 0.5]
        # Compute disordered fraction
        dis_content = len(dis_res) / len(self.values)
        return dis_content

    def continuous_disorder(self):
        """
        Store continuous disordered (CD) regions
        """
        # Filter residues with score equal or greater than 0.5
        # and store its position index
        dis_res = [index for index, res in enumerate(self.values)
                   if float(res) >= 0.5]
        # Count of residues in regions of 2 or more amino acid-long
        cont_dis = 0
        # Counter to store partial results of each continuous region
        c = 0
        # Iterate over disordered residues list
        for i, j in zip(dis_res, dis_res[1:]):
            # Check if residues are consecutive
            if j - i == 1:
                # Update counter
                c += 1
            # Not consecutive
            else:
                # Compute regions with 2 or more residues
                if c > 1:
                    # Add last residue of the interval
                    c += 1
                    # Update disorder count
                    cont_dis += c
                # Reset counter for the next interval
                c = 1
        # Evaluate the last position, as zip only works until
        # the penultimate position
        # if dis_res[-1] - dis_res[-2] == 1:
        #    cont_dis += 1
        CD_fraction = cont_dis / len(self.values)
        return CD_fraction

    def longest_continuous_disorder(self):
        """
        Return longest continuous disorder (CDl) lenght
        and the CDl percentage of protein lenght.
        """
        # Filter residues with score equal or greater than 0.5
        # and store its position index
        dis_res = [index for index, res in enumerate(self.values)
                   if float(res) >= 0.5]
        # Initialize longest CD region
        CDl = 0
        # Counter to store partial results of each continuous region
        c = 0
        # Iterate over disordered residues list
        for i, j in zip(dis_res, dis_res[1:]):
            # Check if residues are consecutive
            if j - i == 1:
                # Update counter
                c += 1
            # Not consecutive
            else:
                # Add last residue of the interval
                c += 1
                # Update CDl
                if c > CDl:
                    CDl = c
                # Reset counter for the next interval
                c = 1
        CDl_fraction = CDl / len(self.values)
        return CDl_fraction


def parse_iupred2a(output):
    """
    Parse IUPRED2A output into a dictionary with
    IUPRED2A and ANCHOR keys.
    """

    predictions = {}
    with open(output, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                identifier = line.strip()[1:]
                predictions[identifier] = {
                                    "iupred2a": [],
                                    "anchor": []
                        }
            if line[0] not in [">", "#"]:
                predictions[identifier]["iupred2a"].append(line.split()[2])
                predictions[identifier]["anchor"].append(line.split()[3])
    return predictions

#    print(results["lcl|NC_038371.1_prot_ORF_96972:96613_predicted"]["iupred2a"])
#    print(len(results["lcl|NC_038371.1_prot_ORF_96972:96613_predicted"]["iupred2a"]))
#    tmp = DisorderFeatures(results["lcl|NC_038371.1_prot_ORF_96972:96613_predicted"]["iupred2a"])
#    print(tmp.disorder_content())
#    print(tmp.continuous_disorder())
#    print(tmp.longest_continuous_disorder())


def parse_disembl(output):
    """
    Parse DisEMBL output into a dictionary with
    HOTLOOPS and REMARK465 keys.
    """

    predictions = {}
    with open(output, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                identifier = line.strip()[1:]
                predictions[identifier] = {
                                    "hotloops": [],
                                    "remark465": []
                        }
            if line[0] not in [">", "#"]:
                predictions[identifier]["hotloops"].append(line.split()[3])
                predictions[identifier]["remark465"].append(line.split()[4])
    return predictions


def compute_features(predictions, predictor):
    """
    """
    if predictor == "iupred2a":
        features = {
                "Protein": [],
                "Disorder_IUPRED2A": [],
                "CD_IUPRED2A": [],
                "CDl_IUPRED2A": [],
                "Disorder_ANCHOR": [],
                "CD_ANCHOR": [],
                "CDl_ANCHOR": [],
                }
        for key in predictions:
            features["Protein"].append(key)

            handle = DisorderFeatures(predictions[key]["iupred2a"])
            features["Disorder_IUPRED2A"].append(handle.disorder_content())
            features["CD_IUPRED2A"].append(handle.continuous_disorder())
            features["CDl_IUPRED2A"].append(handle.longest_continuous_disorder())

            handle = DisorderFeatures(predictions[key]["anchor"])
            features["Disorder_ANCHOR"].append(handle.disorder_content())
            features["CD_ANCHOR"].append(handle.continuous_disorder())
            features["CDl_ANCHOR"].append(handle.longest_continuous_disorder())

    if predictor == "disembl:":
        features = {
                "Protein": [],
                "Disorder_HOTLOOPS": [],
                "CD_HOTLOOPS": [],
                "CDl_HOTLOOPS": [],
                "Disorder_REMARK465": [],
                "CD_REMARK465": [],
                "CDl_REMARK465": [],
                }
        for key in predictions:
            features["Protein"].append(key)

            handle = DisorderFeatures(predictions[key]["hotloops"])
            features["Disorder_HOTLOOPS"].append(handle.disorder_content())
            features["CD_HOTLOOPS"].append(handle.continuous_disorder())
            features["CDl_HOTLOOPS"].append(handle.longest_continuous_disorder())

            handle = DisorderFeatures(predictions[key]["remark465"])
            features["Disorder_REMARK465"].append(handle.disorder_content())
            features["CD_REMARK465"].append(handle.continuous_disorder())
            features["CDl_REMARK465"].append(handle.longest_continuous_disorder())
    df = pd.DataFrame(features)
    df.set_index("Protein")
    df.to_csv("tmp", index=False)


def main():
    """
    Command line argument parser
    """
    parser = argparse.ArgumentParser(
            usage="python3 disorder_features.py <output> <program>",
            description="""
            Compute disorder features from IUPRED2A or DisEMBL predictions
            """,
            epilog=""
    )
    parser.add_argument("output", help="Output file")
    parser.add_argument("program",
                        type=str,
                        choices=["iupred2a", "disembl"],
                        help="Program used for disorder prediction")
    args = parser.parse_args()
    if args.program == "iupred2a":
        predictions = parse_iupred2a(args.output)
    elif args.program == "disembl":
        predictions = parse_disembl(args.output)
    compute_features(predictions, args.program)


if __name__ == "__main__":
    main()
