class DisorderFeatures():
    pass


def parse_iupred2a_output(output):
    """
    Parse IUPRED2A output into a dictionary with
    IUPRED2A and ANCHOR keys.
    """

    results = {}
    with open(output, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                identifier = line.strip()[1:]
                results[identifier] = {
                                    "IUPRED2A": [],
                                    "ANCHOR": []
                        }
            if line[0] not in [">", "#"]:
                results[identifier]["IUPRED2A"].append(line.split()[2])
                results[identifier]["ANCHOR"].append(line.split()[3])
    print(results)
    print(len(results))


def parse_disembl_output(output):
    """
    Parse DisEMBL output into a dictionary with
    COILS and HOTLOOPS keys.
    """

    results = {}
    with open(output, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                identifier = line.strip()[1:]
                results[identifier] = {
                                    "COILS": [],
                                    "HOTLOOPS": []
                        }
            if line[0] not in [">", "#"]:
                results[identifier]["COILS"].append(line.split()[2])
                results[identifier]["HOTLOOPS"].append(line.split()[3])
    print(results)
    print(len(results))


parse_iupred2a_output("CDS-space_iupred2a.txt")
parse_disembl_output("CDS-space_disembl.txt")
