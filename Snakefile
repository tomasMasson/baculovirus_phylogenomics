DATA_DIR = config['data_dir']
ORTHOF_DIR = config['orthof_dir']
MIN_NUMBER = config['min_number']

rule all:
    input:
        "results/hmm_report.csv"

include: 'rules/create_orthogroups.smk'
