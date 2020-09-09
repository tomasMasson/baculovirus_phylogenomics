import glob

SPECIES = [file.split('/')[-1].split('.')[0] for file in glob.glob('data/genomes/*')]
PROTEOMES = [file.split('/')[-1].split('.')[0] for file in glob.glob('data/proteomes/*')]

#SPECIES = ['NC_001623', 'NC_008348']

#rule all:
#  input:
#    "results/CDS_annotation.tbl"

rule CDS_prediction:
  input:
    "data/genomes/{specie}.fa"
  output:
    temp("data/proteomes/{specie}.CDS.fasta")
  shell:
    """
    bin/ORFfinder -in {input} -out {output} -s 0 -ml 150 -n t -c t
    """

rule aggregate_predicted_CDS:
  input:
    expand("data/proteomes/{species}.CDS.fasta", species=SPECIES)
  output:
    temp("results/raw_predicted_CDS.faa")
  shell:
    """
    cat {input} > {output}
    """

rule rename_predicted_CDS:
  input:
    "results/raw_predicted_CDS.faa"
  output:
    temp("results/predicted_CDS.faa")
  shell:
    """
    scripts/rename_predicted_CDS.py {input} > {output}
    """

rule aggregate_annotated_CDS:
  input:
    expand("data/proteomes/{species}.fasta", species=SPECIES)
  output:
    temp("results/annotated_CDS.faa")
  shell:
    """
    cat {input} > {output}
    """

rule merge_CDS:
  input:
    "results/predicted_CDS.faa",
    "results/annotated_CDS.faa"
  output:
    "results/baculovirus_CDS.faa"
  shell:
    """
    scripts/annotate_CDS.py {input} {output}
    """

rule pfam_annotation:
  input:
    seq="results/baculovirus_CDS.faa",
    db="data/Pfam/Pfam-A.hmm",
  output:
    "results/CDS_annotation.tbl"
  threads:
      4
  shell:
    """
    hmmscan --tblout {output} -E 0.001 --domE 0.001 --cpu {threads} {input.db} {input.seq}
    """
