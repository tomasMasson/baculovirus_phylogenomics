#!/usr/bin/env python3

import argparse
import pandas as pd


def longest_CD(values):
    """
    Return the sequence range for the longest continuous
    disorder (CDl) subsequence.
    """
    # Filter residues with score equal or greater than 0.5
    # and store its position index
    dis_res = [index for index, res in enumerate(values)
               if float(res) >= 0.5]
    # Initialize longest CD region
    CDl = []
    # Counter to store partial results of each continuous region
    c = []
    # Iterate over disordered residues list
    for i, j in zip(dis_res, dis_res[1:]):
        # Check if residues are consecutive
        if j - i == 1:
            # Update counter
            c.append(i)
        # Not consecutive
        else:
            # Add last residue of the interval
            c.append(i)
            # Update CDl
            if len(c) > len(CDl):
                CDl = c
            # Reset counter for the next interval
            c = []
    return CDl


class DisorderFeatures:
    def __init__(self, iupred, anchor):
        """
        Initializes with a list of IUPred2A and ANCHOR
        predicted values.
        """
        self.iupred = iupred
        self.anchor = anchor

    def disorder_content(self):
        """
        Filter disordered residues with a predicted
        value greater than 0.5 and compute average
        disorder content
        """
        # Filter residues with score equal or greater than 0.5
        dis_res = [res for res in self.iupred if float(res) >= 0.5]
        # Compute disordered fraction
        dis_content = len(dis_res) / len(self.iupred) * 100
        return dis_content

    def CDl_fraction(self):
        """
        Computes the protein percentage occupied by CDl.
        """
        CDl = longest_CD(self.iupred)
        CDl_fraction = len(CDl) / len(self.iupred) * 100
        return CDl_fraction

    def CDl_lenght(self):
        """
        Returns CDl amino acid lenght.
        """
        return len(longest_CD(self.iupred))

    def CDl_position(self):
        """
        Return CDl centroid position in protein sequence.
        """
        CDl = longest_CD(self.iupred)
        if CDl:
            return sum(CDl) / len(CDl)
        else:
            return 0

    def anchor_position(self):
        """
        Return ANCHOR centroid position in protein sequence.
        """
        anchor = longest_CD(self.anchor)
        if anchor:
            return sum(anchor) / len(anchor)
        else:
            return 0


def parse_iupred2a(output):
    """
    Parse IUPRED2A output into a dictionary with
    IUPRED2A and ANCHOR predictions.
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


# def parse_disembl(output):
#    """
#    Parse DisEMBL output into a dictionary with
#    HOTLOOPS and REMARK465 keys.
#    """
#
#    predictions = {}
#    with open(output, "r") as fh:
#        for line in fh:
#            if line.startswith(">"):
#                identifier = line.strip()[1:]
#                predictions[identifier] = {
#                                    "hotloops": [],
#                                    "remark465": []
#                        }
#            if line[0] not in [">", "#"]:
#                predictions[identifier]["hotloops"].append(line.split()[3])
#                predictions[identifier]["remark465"].append(line.split()[2])
#    return predictions


def compute_features(predictions):
    """
    Generates a CSV file with disorder features for each
    protein present in the input file.
    """
#    if predictor == "iupred2a":
    features = {
            "Protein": [],
            "Disorder_content": [],
            "LCPL": [],
            "CDl": [],
            "CDl_position": [],
            "ANCHOR_position": [],
            }
    for key in predictions:
        protein_id = key

        # Create DisorderFeatures instance
        iupred = predictions[key]["iupred2a"]
        anchor = predictions[key]["anchor"]
        handle = DisorderFeatures(iupred, anchor)

        # Compute features
        dis_content = handle.disorder_content()
        LCPL = handle.CDl_fraction()
        CDl = handle.CDl_lenght()
        CDl_position = handle.CDl_position()
        ANCHOR_position = handle.anchor_position()

        # Store features
        features["Protein"].append(protein_id)
        features["Disorder_content"].append(dis_content)
        features["LCPL"].append(LCPL)
        features["CDl"].append(CDl)
        features["CDl_position"].append(CDl_position)
        features["ANCHOR_position"].append(ANCHOR_position)
    df = pd.DataFrame(features)
    df.set_index("Protein")
    df.to_csv("tmp", index=False)

#    if predictor == "disembl":
#        features = {
#                "Protein": [],
#                "Disorder_HOTLOOPS": [],
#                "CD_HOTLOOPS": [],
#                "CDl_HOTLOOPS": [],
#                "Disorder_REMARK465": [],
#                "CD_REMARK465": [],
#                "CDl_REMARK465": [],
#                }
#        for key in predictions:
#            features["Protein"].append(key)
#
#            handle = DisorderFeatures(predictions[key]["hotloops"])
#            features["Disorder_HOTLOOPS"].append(handle.disorder_content())
#            features["CD_HOTLOOPS"].append(handle.continuous_disorder())
#            features["CDl_HOTLOOPS"].append(handle.longest_continuous_disorder())
#
#            handle = DisorderFeatures(predictions[key]["remark465"])
#            features["Disorder_REMARK465"].append(handle.disorder_content())
#            features["CD_REMARK465"].append(handle.continuous_disorder())
#            features["CDl_REMARK465"].append(handle.longest_continuous_disorder())


def main():
    """
    Command line argument parser
    """
    parser = argparse.ArgumentParser(
            usage="python3 disorder_features.py <output>",
            description="""
            Compute disorder features from IUPRED2A predictions
            """,
            epilog=""
    )
    parser.add_argument("output", help="Output file")
#    parser.add_argument("program",
#                        type=str,
#                        choices=["iupred2a", "disembl"],
#                        help="Program used for disorder prediction")
    args = parser.parse_args()
#    if args.program == "iupred2a":
#        predictions = parse_iupred2a(args.output)
#    elif args.program == "disembl":
#        predictions = parse_disembl(args.output)
    predictions = parse_iupred2a(args.output)
    compute_features(predictions)


if __name__ == "__main__":
    main()
