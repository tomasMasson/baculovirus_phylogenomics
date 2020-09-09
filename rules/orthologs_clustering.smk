from datetime import datetime

date = datetime.now().strftime('%h') + datetime.now().strftime('%d')
# Set path to Orthofinder orthogroup sequences data
ORTHOF_DIR = os.path.join(f'OrthoFinder/Results_{date}', 'Orthogroup_Sequences')

#rule all:
#  input:
#    "results/gene_clusters.csv"

rule split_ORFomes:
    input:
      "results/baculovirus_CDS.faa"
    output:
      directory("results/ORFomes/")
    shell:
        """
        scripts/split_ORFomes.py {input} &&\
        mv NC_*.faa {output}
        """

rule run_orthofinder:
    input:
        "results/ORFomes/"
    output:
        "results/Orthogroup_Sequences/OG0000000.fa"
    threads:
        4
    shell:
        """
        orthofinder -I 1.4 -t 4 -f {input} -og -o ./OrthoFinder &&\
        cp -r OrthoFinder {ORTHOF_DIR} results &&\
        rm -rf OrthoFinder/
        """

rule create_hmm_database:
    input:
        "results/Orthogroup_Sequences/OG0000000.fa"
    output:
        "results/hmm_db"
    params:
        "results/Orthogroup_Sequences",
        min_numer=6,
        name='hmm_db'
    shell:
        "python3 scripts/create_hmm_database.py {params} > {output}"

rule hmmer_search:
    input:
        seqdb="results/baculovirus_CDS.faa",
        hmmfile="results/hmm_db"
    output:
        "results/hmm_search.tab"
    shell:
        """
        hmmscan --tblout {output} -E 0.001 --domE 0.001 {input.hmmfile} {input.seqdb}
        """

rule hmmer_report:
    input:
        "results/hmm_search.tab",
        "data/protein_names.csv"
    output:
        "results/gene_clusters.csv"
    shell:
        """
        python3 scripts/parse_hmmer-tab.py {input} > {output}
        """
