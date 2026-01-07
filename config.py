# Configuration for dehalogenases investigation

# PDB structure to analyze
PDB_ID = "3BJX"

# BLAST search parameters
BLAST_DATABASE = "swissprot"  # UniProt Swiss-Prot database
E_VALUE_THRESHOLD = 0.001
MAX_HITS = 100

# MSA parameters
MSA_TOOL = "clustalo"  # Clustal Omega
MSA_OUTPUT_FORMAT = "fasta"

# Output directories
DATA_DIR = "data"
RESULTS_DIR = "results"

# File paths
PDB_SEQUENCE_FILE = f"{DATA_DIR}/3BJX_sequence.fasta"
BLAST_RESULTS_FILE = f"{DATA_DIR}/blast_results.xml"
HOMOLOGS_FILE = f"{DATA_DIR}/homologs.fasta"
MSA_FILE = f"{RESULTS_DIR}/alignment.fasta"
MSA_VISUALIZATION_FILE = f"{RESULTS_DIR}/msa_visualization.html"
