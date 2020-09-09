rule all:
    input:
      "results/baculovirus_CDS.faa",
      "results/gene_clusters.csv"

include: "rules/CDS-space_prediction.smk"
include: "rules/orthologs_clustering.smk"
