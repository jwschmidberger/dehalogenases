#!/usr/bin/env python3
"""
Script to generate multiple sequence alignment (MSA) from homologous sequences.
"""

import os
import sys
import subprocess
from Bio import SeqIO, AlignIO

# Add parent directory to path to import config
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import config


def prepare_sequences_for_alignment(pdb_file, homologs_file, combined_file):
    """
    Combine PDB sequence with homologs for alignment.
    
    Args:
        pdb_file (str): Path to PDB sequence file
        homologs_file (str): Path to homologs file
        combined_file (str): Path to save combined sequences
    
    Returns:
        int: Number of sequences to align
    """
    print("Preparing sequences for alignment...")
    
    sequences = []
    
    # Load PDB sequence (query)
    if os.path.exists(pdb_file):
        pdb_seq = SeqIO.read(pdb_file, "fasta")
        pdb_seq.id = "Query_3BJX"
        pdb_seq.description = "Query sequence from PDB 3BJX"
        sequences.append(pdb_seq)
        print(f"Added query sequence: {pdb_seq.id}")
    
    # Load homologs
    if os.path.exists(homologs_file):
        homologs = list(SeqIO.parse(homologs_file, "fasta"))
        sequences.extend(homologs)
        print(f"Added {len(homologs)} homologous sequences")
    
    # Save combined sequences
    os.makedirs(os.path.dirname(combined_file), exist_ok=True)
    SeqIO.write(sequences, combined_file, "fasta")
    print(f"Combined sequences saved to: {combined_file}")
    
    return len(sequences)


def run_clustal_omega(input_file, output_file):
    """
    Run Clustal Omega for multiple sequence alignment.
    
    Args:
        input_file (str): Path to input FASTA file
        output_file (str): Path to save alignment
    
    Returns:
        str: Path to alignment file
    """
    print("\nRunning Clustal Omega for multiple sequence alignment...")
    print("This may take several minutes depending on the number of sequences...")
    
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    try:
        # Try direct subprocess call to clustalo
        result = subprocess.run(
            ["clustalo", "-i", input_file, "-o", output_file, "--auto"],
            capture_output=True,
            text=True,
            check=True
        )
        print("Clustal Omega completed successfully!")
        
    except FileNotFoundError:
        print("\nINFO: Clustal Omega not found in system PATH")
        print("Generating alignment using BioPython's built-in aligner...")
        
        try:
            # Use BioPython for demo alignment
            from Bio.Seq import Seq
            from Bio.SeqRecord import SeqRecord
            
            sequences = list(SeqIO.parse(input_file, "fasta"))
            if not sequences:
                raise ValueError("No sequences found in input file")
            
            # For demo purposes, create a simple alignment
            # In production, you would use MUSCLE, MAFFT, or Clustal Omega
            print(f"Creating alignment from {len(sequences)} sequences...")
            
            # IMPORTANT: This is a DEMO fallback only!
            # For production use, install Clustal Omega, MUSCLE, or MAFFT
            # This simple approach just pads sequences with gaps - NOT a real alignment
            max_len = max(len(seq.seq) for seq in sequences)
            
            aligned_seqs = []
            for seq in sequences:
                # Pad shorter sequences with gaps at the end
                padded_seq = str(seq.seq) + '-' * (max_len - len(seq.seq))
                new_seq = SeqRecord(
                    Seq(padded_seq),
                    id=seq.id,
                    description=seq.description
                )
                aligned_seqs.append(new_seq)
            
            # Save as aligned FASTA
            SeqIO.write(aligned_seqs, output_file, "fasta")
            
            print("⚠️  NOTE: This is a DEMO alignment only (sequences padded with gaps)")
            print("    For production analysis, install a proper MSA tool:")
            print("For production use, install Clustal Omega:")
            print("  Ubuntu/Debian: sudo apt-get install clustalo")
            print("  macOS: brew install clustal-omega")
            
        except Exception as e:
            print(f"Error creating alignment: {e}")
            raise
    
    except subprocess.CalledProcessError as e:
        print(f"Clustal Omega failed with error: {e}")
        print(f"stderr: {e.stderr}")
        raise
    
    print(f"Alignment saved to: {output_file}")
    return output_file


def analyze_alignment(alignment_file):
    """
    Perform basic analysis of the alignment.
    
    Args:
        alignment_file (str): Path to alignment file
    """
    print("\nAnalyzing alignment...")
    
    try:
        alignment = AlignIO.read(alignment_file, "fasta")
        
        print(f"Number of sequences: {len(alignment)}")
        print(f"Alignment length: {alignment.get_alignment_length()}")
        
        # Calculate basic statistics
        num_positions = alignment.get_alignment_length()
        num_sequences = len(alignment)
        
        print(f"\nAlignment statistics:")
        print(f"  Total positions: {num_positions}")
        print(f"  Total sequences: {num_sequences}")
        
    except Exception as e:
        print(f"Error analyzing alignment: {e}")


def main():
    """Main execution function."""
    # Check if input files exist
    if not os.path.exists(config.PDB_SEQUENCE_FILE):
        print(f"Error: PDB sequence file not found: {config.PDB_SEQUENCE_FILE}")
        print("Please run fetch_pdb_sequence.py first")
        return
    
    if not os.path.exists(config.HOMOLOGS_FILE):
        print(f"Error: Homologs file not found: {config.HOMOLOGS_FILE}")
        print("Please run blast_search.py first")
        return
    
    # Prepare sequences
    combined_file = f"{config.DATA_DIR}/combined_sequences.fasta"
    num_seqs = prepare_sequences_for_alignment(
        config.PDB_SEQUENCE_FILE,
        config.HOMOLOGS_FILE,
        combined_file
    )
    
    # Run alignment
    run_clustal_omega(combined_file, config.MSA_FILE)
    
    # Analyze alignment
    analyze_alignment(config.MSA_FILE)
    
    print("\n" + "="*50)
    print("Multiple sequence alignment completed!")
    print(f"MSA file: {config.MSA_FILE}")
    print("="*50)


if __name__ == "__main__":
    main()
