DATA_DIR = config['data_dir']
OG_DIR = config['og_dir']

rule all:
    input:
        "results/hmm_report.csv"

#rule run_broccoli:
#    input:
#        DATA_DIR
#    output:
#        directory("results/Broccoli/")
#    shell:
#        """
#        python3 ~/Broccoli/broccoli.py -dir {input} -path_diamond ~/diamond &&\
#        mv dir_step* {output}
#        """
#
rule run_orthofinder:
    input:
        DATA_DIR
    output:
        directory("results/Orthofinder/")
    shell:
        """
        orthofinder -f {input} -og -o ./Orthofinder &&\
        mv Orthofinder results
        """

rule build_hmm_database:
    input:
        rules.run_orthofinder.output[0]
    output:
        "results/hmmdb"
    params:
        os.path.join(rules.run_orthofinder.output[0], OG_DIR)
    shell:
        "scripts/generate_hmm_database.py {params} > {output}"

#rule aggregate_protein_sequences:
#    input:
#    output:
#        "results/proteins_db.faa"
#    params:
#        DATA_DIR
#    shell:
#        """
#        cat {params}/* > {output}
#        """

rule hmmer_search:
    input:
        seqdb="data/putative_proteins_baculovirus.faa",
        hmmfile="results/hmmdb"
    output:
        "results/hmm_search.out"
    shell:
        """
        hmmscan -o {output} -E 0.001 --domE 0.001 {input.hmmfile} {input.seqdb}
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
