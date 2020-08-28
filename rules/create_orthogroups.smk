DATA_DIR = config['data_dir']
ORTHOF_DIR = config['orthof_dir']
MIN_NUMBER = config['min_number']

rule run_orthofinder:
    input:
        DATA_DIR
    output:
        "results/Orthogroup_Sequences/OG0000000.fa"
    threads:
        3
    shell:
        """
        orthofinder -I 1.4 -t 3 -f {input} -og -o ./OrthoFinder &&\
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
        MIN_NUMBER,
        name='hmm_db'
    shell:
        "python3 scripts/create_hmm_database.py {params} > {output}"

rule aggregate_protein_sequences:
    input:
    output:
        temp("results/proteins_db.faa")
    params:
        DATA_DIR
    shell:
        """
        cat {params}/* > {output}
        """

rule hmmer_search:
    input:
        seqdb="results/proteins_db.faa",
        hmmfile="results/hmm_db"
    output:
        "results/hmm_search.tab"
    shell:
        """
        hmmscan --domtblout {output} -E 0.001 --domE 0.001 {input.hmmfile} {input.seqdb}
        """

rule hmmer_report:
    input:
        "results/hmm_search.tab"
    output:
        "results/gene_clusters.csv"
    shell:
        """
        python3 scripts/parse_hmmer-tab.py {input} > {output}
        """
