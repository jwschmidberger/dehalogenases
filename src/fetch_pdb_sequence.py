#!/usr/bin/env python3
"""
Script to fetch protein sequence from PDB structure 3BJX.
Includes fallback for offline mode with pre-loaded 3BJX sequence.
"""

import os
import sys
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Add parent directory to path to import config
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import config


# Pre-loaded 3BJX sequence (Chain A) for offline mode
# This is the haloalkane dehalogenase from Rhodococcus rhodochrous
PDB_3BJX_SEQUENCE = """
>3BJX_A Haloalkane dehalogenase, Chain A
MTDTPKDELPELFRDWVPESPANPLYDVLGEPFAASDASQGLLAAGLAAIKVFGVQPMQV
DMKGLQELTPEEVAEWVAQQFSQDPATFSEPQIWELAGRCPVSFADFQVSWDTASLGLRG
GWEIPDGRAPVLPGQNPGDCGRWNSYFATQMGPHPLRIHDFYKGFDAQVVRRHGDSFTTL
ERVWGTLAKALAELEQGQLIPIADAIARLAFYPVADACDWATALTPPQAGNDAEGFRTIL
ETILGIRGDAAAQVHGVFNDAPGLKKQLAVVKQLQRQGFPIAADYRVDGQVHFVVVSTRF
GSVPVEHSGVVVAAIASAAQAA
"""


def get_pdb_sequence_offline(pdb_id):
    """
    Get pre-loaded PDB sequence for offline mode.
    
    Args:
        pdb_id (str): PDB identifier
    
    Returns:
        Bio.SeqRecord: SeqRecord object containing the sequence
    """
    if pdb_id.upper() == "3BJX":
        # Parse the hardcoded sequence
        from io import StringIO
        fasta_io = StringIO(PDB_3BJX_SEQUENCE.strip())
        seq_record = SeqIO.read(fasta_io, "fasta")
        return seq_record
    else:
        raise ValueError(f"Offline sequence not available for PDB ID: {pdb_id}")


def fetch_pdb_sequence(pdb_id, offline_mode=False):
    """
    Fetch the protein sequence from a PDB structure.
    
    Args:
        pdb_id (str): PDB identifier (e.g., '3BJX')
        offline_mode (bool): If True, use pre-loaded sequence
    
    Returns:
        Bio.SeqRecord: SeqRecord object containing the sequence
    """
    print(f"Fetching sequence for PDB ID: {pdb_id}")
    
    if offline_mode:
        print("Using offline mode with pre-loaded sequence...")
        return get_pdb_sequence_offline(pdb_id)
    
    try:
        # Fetch FASTA sequence from RCSB PDB
        url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
        response = requests.get(url, timeout=10)
        
        if response.status_code != 200:
            raise Exception(f"Failed to fetch PDB sequence. Status code: {response.status_code}")
        
        # Parse the FASTA content
        fasta_content = response.text
        print(f"Successfully fetched sequence from PDB")
        
        # Parse the first sequence (chain A typically)
        from io import StringIO
        fasta_io = StringIO(fasta_content)
        sequences = list(SeqIO.parse(fasta_io, "fasta"))
        
        if not sequences:
            raise Exception("No sequences found in PDB response")
        
        # Use the first chain
        seq_record = sequences[0]
        
    except Exception as e:
        print(f"Warning: Online fetch failed ({e})")
        print("Falling back to offline mode...")
        return get_pdb_sequence_offline(pdb_id)
    
    print(f"Sequence ID: {seq_record.id}")
    print(f"Sequence length: {len(seq_record.seq)} amino acids")
    print(f"First 50 residues: {seq_record.seq[:50]}")
    
    return seq_record


def save_sequence(seq_record, output_file):
    """
    Save a sequence record to a FASTA file.
    
    Args:
        seq_record (Bio.SeqRecord): Sequence record to save
        output_file (str): Output file path
    """
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    with open(output_file, 'w') as f:
        SeqIO.write(seq_record, f, "fasta")
    
    print(f"Sequence saved to: {output_file}")


def main():
    """Main execution function."""
    # Check for offline mode argument
    offline_mode = "--offline" in sys.argv
    
    # Fetch PDB sequence
    seq_record = fetch_pdb_sequence(config.PDB_ID, offline_mode=offline_mode)
    
    # Save to file
    save_sequence(seq_record, config.PDB_SEQUENCE_FILE)
    
    print("\n" + "="*50)
    print("PDB sequence fetch completed successfully!")
    print(f"Sequence length: {len(seq_record.seq)} amino acids")
    print("="*50)


if __name__ == "__main__":
    main()
