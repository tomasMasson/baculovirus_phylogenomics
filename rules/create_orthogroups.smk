#DATA_DIR = config['data_dir']
#ORTHOF_DIR = config['orthof_dir']
#MIN_NUMBER = config['min_number']
#
#rule all:
#    input:
#        "results/hmm_db",
#        "results/hmm_report.csv"
#
rule run_orthofinder:
    input:
        DATA_DIR
    output:
        "results/Orthogroup_Sequences/OG0000000.fa"
    shell:
        """
        orthofinder -f {input} -og -o ./OrthoFinder &&\
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
        "results/proteins_db.faa"
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
        "results/hmm_search.out"
    shell:
        """
        hmmscan -o {output} -E 0.01 --domE 0.01 {input.hmmfile} {input.seqdb}
        """

rule hmmer_report:
    input:
        "results/hmm_search.out"
    output:
        "results/hmm_report.csv"
    shell:
        """
        python3 scripts/extract_hmmer_hits.py {input} > {output}
        """
