import glob

SPECIES = [file.split('/')[-1].split('.')[0] for file in glob.glob('data/genomes/*')]
PROTEOMES = [file.split('/')[-1].split('.')[0] for file in glob.glob('data/proteomes/*')]

#SPECIES = ['NC_001623', 'NC_008348']

rule all:
  input:
    "results/pfam_annotation.tbl"

rule orf_prediction:
  input:
    "data/genomes/{specie}.fa"
  output:
    temp("data/proteomes/{specie}.orffinder.fasta")
  shell:
    """
    scripts/ORFfinder -in {input} -out {output} -s 0 -ml 150 -n t -c t
    """

rule aggregate_predicted_orfs:
  input:
    expand("data/proteomes/{species}.orffinder.fasta", species=SPECIES)
  output:
    temp("results/baculovirus_predicted_orfs.faa")
  shell:
    """
    cat {input} > {output}
    """

rule aggregate_annotated_orfs:
  input:
    expand("data/proteomes/{species}.fasta", species=SPECIES)
  output:
    temp("results/baculovirus_annotated_orfs.faa")
  shell:
    """
    cat {input} > {output}
    """

rule rename_predicted_orfs:
  input:
    "results/baculovirus_predicted_orfs.faa",
    "results/baculovirus_annotated_orfs.faa"
  output:
    "results/baculovirus_orfs.faa"
  shell:
    "scripts/rename_predicted_orfs.py {input} {output}"

rule pfam_annotation:
  input:
    seq="results/baculovirus_orfs.faa",
    db="data/Pfam/Pfam-A.hmm",
  output:
    "results/pfam_annotation.tbl"
  threads:
      4
  shell:
    """
    hmmscan --tblout {output} -E 0.001 --domE 0.001 --cpu {threads} {input.db} {input.seq}
    """
