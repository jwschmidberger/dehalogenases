#!/usr/bin/env python3
"""
Script to perform BLAST search using the 3BJX sequence to find homologous sequences.
"""

import os
import sys
import time
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

# Add parent directory to path to import config
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import config


def run_blast_search(sequence_file, output_file):
    """
    Run BLAST search using NCBI's online BLAST service.
    
    Args:
        sequence_file (str): Path to input FASTA file
        output_file (str): Path to save BLAST results (XML format)
    
    Returns:
        str: Path to BLAST results file
    """
    print("Loading query sequence...")
    seq_record = SeqIO.read(sequence_file, "fasta")
    sequence = str(seq_record.seq)
    
    print(f"Running BLAST search against {config.BLAST_DATABASE}...")
    print(f"Query length: {len(sequence)} amino acids")
    print(f"E-value threshold: {config.E_VALUE_THRESHOLD}")
    print(f"Max hits: {config.MAX_HITS}")
    print("\nThis may take several minutes...")
    
    # Run BLAST search
    result_handle = NCBIWWW.qblast(
        program="blastp",
        database=config.BLAST_DATABASE,
        sequence=sequence,
        expect=config.E_VALUE_THRESHOLD,
        hitlist_size=config.MAX_HITS
    )
    
    print("BLAST search completed!")
    
    # Save results
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, 'w') as out_handle:
        out_handle.write(result_handle.read())
    
    result_handle.close()
    print(f"BLAST results saved to: {output_file}")
    
    return output_file


def parse_blast_results(blast_file, output_fasta):
    """
    Parse BLAST results and extract homologous sequences.
    
    Args:
        blast_file (str): Path to BLAST XML results
        output_fasta (str): Path to save homologous sequences
    
    Returns:
        list: List of hit accessions
    """
    print("\nParsing BLAST results...")
    
    with open(blast_file) as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        blast_record = next(blast_records)
    
    # Collect accessions from significant hits
    accessions = []
    hits_data = []
    
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < config.E_VALUE_THRESHOLD:
                # Extract accession from hit definition
                hit_def = alignment.hit_def
                hit_id = alignment.hit_id
                accession = hit_id.split('|')[1] if '|' in hit_id else hit_id
                
                hits_data.append({
                    'accession': accession,
                    'description': hit_def,
                    'e_value': hsp.expect,
                    'identity': hsp.identities,
                    'align_length': hsp.align_length,
                    'sequence': hsp.sbjct
                })
                accessions.append(accession)
    
    print(f"Found {len(accessions)} significant hits")
    
    # Save hit information
    if hits_data:
        os.makedirs(os.path.dirname(output_fasta), exist_ok=True)
        
        with open(output_fasta, 'w') as out_handle:
            for i, hit in enumerate(hits_data[:config.MAX_HITS], 1):
                # Write in FASTA format
                header = f">{hit['accession']} {hit['description']}"
                out_handle.write(f"{header}\n")
                # Write sequence in lines of 60 characters
                sequence = hit['sequence'].replace('-', '')  # Remove gaps
                for j in range(0, len(sequence), 60):
                    out_handle.write(f"{sequence[j:j+60]}\n")
        
        print(f"Homologous sequences saved to: {output_fasta}")
    
    return accessions


def main():
    """Main execution function."""
    # Check if sequence file exists
    if not os.path.exists(config.PDB_SEQUENCE_FILE):
        print(f"Error: Sequence file not found: {config.PDB_SEQUENCE_FILE}")
        print("Please run fetch_pdb_sequence.py first")
        return
    
    # Run BLAST search
    blast_file = run_blast_search(config.PDB_SEQUENCE_FILE, config.BLAST_RESULTS_FILE)
    
    # Parse results
    accessions = parse_blast_results(blast_file, config.HOMOLOGS_FILE)
    
    print("\n" + "="*50)
    print("BLAST search completed successfully!")
    print(f"Total homologs identified: {len(accessions)}")
    print("="*50)


if __name__ == "__main__":
    main()
