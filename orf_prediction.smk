import glob

SPECIES = [file.split('/')[-1].split('.')[0] for file in glob.glob('data/genomes/*')]

#SPECIES = ['NC_001623', 'NC_008348']

rule all:
  input:
    expand("data/proteomes/{species}.orfs.fasta", species=SPECIES)

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
    "data/proteomes/{specie}.fasta",
    "data/proteomes/{specie}.orffinder.fasta"
  output:
    temp("data/proteomes/{specie}.orfs.aggregated")
  shell:
    """
    cat {input} > {output}
    """

rule deduplicate_orfs:
  input:
    "data/proteomes/{specie}.orfs.aggregated"
  output:
    "data/proteomes/{specie}.orfs.fasta"
  shell:
    """
    cd-hit -i {input} -o {output} -c 1 &&\
    rm {output}.clstr
    """
