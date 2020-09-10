rule all:
    input:
      "results/CDS-space_annotation.tbl",
      "results/gene_clusters.csv"

include: "rules/CDS-space_prediction.smk"
include: "rules/orthologs_clustering.smk"
