rule all:
    input:
      "results/CDS-space_annotation.tbl",
      "results/CDS-space_disembl.txt",
      "results/CDS-space_iupred2a.txt",
      "results/gene_clusters.csv"

include: "rules/CDS-space_prediction.smk"
include: "rules/orthologs_clustering.smk"
