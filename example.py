#!/usr/bin/env python3
"""
Example script demonstrating the complete dehalogenases investigation workflow.
"""

import os
import sys

print("="*60)
print("DEHALOGENASES INVESTIGATION - EXAMPLE WORKFLOW")
print("="*60)
print()

# Step 1: Fetch PDB sequence
print("Step 1: Fetching PDB sequence for 3BJX...")
print("-" * 60)
from src import fetch_pdb_sequence
fetch_pdb_sequence.main()
print()

# Step 2: Check if we have homologous sequences (demo data provided)
print("Step 2: Checking homologous sequences...")
print("-" * 60)
if os.path.exists("data/homologs.fasta"):
    print("✓ Demo homologous sequences available")
    print("  (In production, run blast_search.py to fetch real BLAST results)")
    
    # Count sequences
    from Bio import SeqIO
    homologs = list(SeqIO.parse("data/homologs.fasta", "fasta"))
    print(f"  Found {len(homologs)} homologous dehalogenase sequences")
else:
    print("⚠ No homologous sequences found")
    print("  Run: python src/blast_search.py")
print()

# Step 3: Generate MSA
print("Step 3: Generating Multiple Sequence Alignment...")
print("-" * 60)
from src import generate_msa
generate_msa.main()
print()

# Step 4: Visualize MSA
print("Step 4: Creating visualizations...")
print("-" * 60)
from src import visualize_msa
visualize_msa.main()
print()

# Summary
print("="*60)
print("WORKFLOW COMPLETE!")
print("="*60)
print()
print("Generated files:")
print("  📄 data/3BJX_sequence.fasta - Query sequence from PDB")
print("  📄 data/homologs.fasta - Homologous sequences")
print("  📄 results/alignment.fasta - Multiple sequence alignment")
print("  📊 results/msa_heatmap.png - Alignment heatmap")
print("  📊 results/msa_heatmap_conservation.png - Conservation plot")
print("  📊 results/msa_heatmap_identity.png - Identity matrix")
print()
print("Next steps:")
print("  1. View the alignment file: results/alignment.fasta")
print("  2. Examine the visualizations in the results/ directory")
print("  3. For real BLAST data, run: python src/blast_search.py")
print("     (Note: This requires internet access and takes 5-10 minutes)")
print()
