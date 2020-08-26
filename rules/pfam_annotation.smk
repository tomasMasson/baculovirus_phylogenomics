import glob

SPECIES = [file.split('/')[-1].split('.')[0] for file in glob.glob('data/genomes/*')]

#SPECIES = ['NC_001623', 'NC_008348']

rule all:
  input:
    "pfam_annotation.tbl"
#    expand("data/proteomes/{species}.orfs.fasta", species=SPECIES)

rule orf_prediction:
  input:
    "data/genomes/{specie}.fa"
  output:
    temp("data/proteomes/{specie}.orffinder.fasta")
  shell:
    """
    scripts/ORFfinder -in {input} -out {output} -s 0 -ml 150 -n t -c t
    """

rule aggregate_orfs:
  input:
    expand("data/proteomes/{species}.orffinder.fasta", species=SPECIES)
  output:
    temp("baculovirus_orfs.faa")
  shell:
    """
    cat {input} > {output}
    """

rule pfam_annotation:
  input:
    seq="baculovirus_orfs.faa",
    db="data/Pfam-A.hmm",
  output:
    "pfam_annotation.tbl"
  threads:
      3
  shell:
    """
    hmmscan --tblout {output} -E 0.001 --domE 0.001 --cpu {threads} {input.db} {input.seq}
    """
