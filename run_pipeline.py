#!/usr/bin/env python3
"""
Main pipeline script for dehalogenases investigation.
Runs the complete workflow: fetch PDB sequence, BLAST search, MSA generation, and visualization.
"""

import sys
import os

# Add src directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from src import fetch_pdb_sequence, blast_search, generate_msa, visualize_msa


def print_banner(text):
    """Print a formatted banner."""
    print("\n" + "="*60)
    print(f"  {text}")
    print("="*60 + "\n")


def run_pipeline(skip_blast=False):
    """
    Run the complete dehalogenases investigation pipeline.
    
    Args:
        skip_blast (bool): If True, skip BLAST search (useful for testing)
    """
    print_banner("DEHALOGENASES INVESTIGATION PIPELINE")
    
    try:
        # Step 1: Fetch PDB sequence
        print_banner("STEP 1: Fetching PDB Sequence (3BJX)")
        fetch_pdb_sequence.main()
        
        # Step 2: BLAST search
        if not skip_blast:
            print_banner("STEP 2: Running BLAST Search")
            print("NOTE: This step may take 5-10 minutes...")
            blast_search.main()
        else:
            print_banner("STEP 2: BLAST Search (SKIPPED)")
            print("Using existing BLAST results or demo data")
        
        # Step 3: Generate MSA
        print_banner("STEP 3: Generating Multiple Sequence Alignment")
        generate_msa.main()
        
        # Step 4: Visualize MSA
        print_banner("STEP 4: Creating MSA Visualizations")
        visualize_msa.main()
        
        # Success message
        print_banner("PIPELINE COMPLETED SUCCESSFULLY!")
        print("Results are available in the 'results/' directory:")
        print("  - alignment.fasta: Multiple sequence alignment")
        print("  - msa_heatmap.png: Heatmap visualization")
        print("  - msa_heatmap_conservation.png: Conservation plot")
        print("  - msa_heatmap_identity.png: Sequence identity matrix")
        
    except Exception as e:
        print(f"\n❌ ERROR: Pipeline failed with error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    # Check for command line arguments
    skip_blast = "--skip-blast" in sys.argv
    
    if skip_blast:
        print("⚠️  Running pipeline with BLAST search skipped")
    
    run_pipeline(skip_blast=skip_blast)
