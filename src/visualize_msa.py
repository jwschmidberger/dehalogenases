#!/usr/bin/env python3
"""
Script to visualize multiple sequence alignment (MSA).
"""

import os
import sys
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for headless environments
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap
import numpy as np
from Bio import AlignIO

# Add parent directory to path to import config
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import config


def create_color_scheme():
    """
    Create color scheme for amino acids.
    
    Returns:
        dict: Mapping of amino acids to colors
    """
    # Standard amino acid coloring scheme based on properties
    color_scheme = {
        # Hydrophobic (orange/yellow)
        'A': '#FFA500', 'V': '#FFA500', 'L': '#FFA500', 'I': '#FFA500', 
        'M': '#FFA500', 'F': '#FFA500', 'W': '#FFA500', 'P': '#FFA500',
        
        # Polar (green)
        'S': '#00FF00', 'T': '#00FF00', 'N': '#00FF00', 'Q': '#00FF00', 
        'C': '#00FF00', 'G': '#00FF00', 'Y': '#00FF00',
        
        # Positively charged (blue)
        'K': '#0000FF', 'R': '#0000FF', 'H': '#0000FF',
        
        # Negatively charged (red)
        'D': '#FF0000', 'E': '#FF0000',
        
        # Gap (white)
        '-': '#FFFFFF',
        
        # Unknown (gray)
        'X': '#808080', 'B': '#808080', 'Z': '#808080'
    }
    
    return color_scheme


def visualize_msa_heatmap(alignment_file, output_file, max_sequences=50, max_positions=500):
    """
    Create a heatmap visualization of the MSA.
    
    Args:
        alignment_file (str): Path to alignment file
        output_file (str): Path to save visualization
        max_sequences (int): Maximum number of sequences to display
        max_positions (int): Maximum number of positions to display
    """
    print("Creating MSA heatmap visualization...")
    
    # Read alignment
    alignment = AlignIO.read(alignment_file, "fasta")
    
    num_sequences = len(alignment)
    alignment_length = alignment.get_alignment_length()
    
    print(f"Alignment has {num_sequences} sequences and {alignment_length} positions")
    
    # Limit display size for readability
    display_seqs = min(num_sequences, max_sequences)
    display_pos = min(alignment_length, max_positions)
    
    print(f"Displaying first {display_seqs} sequences and {display_pos} positions")
    
    # Create color scheme
    color_scheme = create_color_scheme()
    
    # Convert alignment to numeric matrix
    aa_list = list(color_scheme.keys())
    aa_to_num = {aa: i for i, aa in enumerate(aa_list)}
    
    matrix = np.zeros((display_seqs, display_pos))
    for i, record in enumerate(alignment[:display_seqs]):
        for j, aa in enumerate(str(record.seq)[:display_pos]):
            aa_upper = aa.upper()
            matrix[i, j] = aa_to_num.get(aa_upper, aa_to_num['X'])
    
    # Create figure
    fig, ax = plt.subplots(figsize=(20, 10))
    
    # Create colormap
    colors = [color_scheme[aa] for aa in aa_list]
    cmap = ListedColormap(colors)
    
    # Plot heatmap
    im = ax.imshow(matrix, cmap=cmap, aspect='auto', interpolation='nearest')
    
    # Set labels
    ax.set_xlabel('Alignment Position', fontsize=12)
    ax.set_ylabel('Sequence', fontsize=12)
    ax.set_title('Multiple Sequence Alignment Heatmap', fontsize=14, fontweight='bold')
    
    # Set y-axis labels (sequence names)
    seq_names = [record.id[:30] for record in alignment[:display_seqs]]  # Truncate long names
    ax.set_yticks(range(display_seqs))
    ax.set_yticklabels(seq_names, fontsize=6)
    
    # Add legend
    legend_elements = [
        mpatches.Patch(color='#FFA500', label='Hydrophobic'),
        mpatches.Patch(color='#00FF00', label='Polar'),
        mpatches.Patch(color='#0000FF', label='Positive'),
        mpatches.Patch(color='#FF0000', label='Negative'),
        mpatches.Patch(color='#FFFFFF', label='Gap'),
        mpatches.Patch(color='#808080', label='Unknown')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10)
    
    plt.tight_layout()
    
    # Save figure
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Heatmap saved to: {output_file}")
    
    plt.close()


def create_conservation_plot(alignment_file, output_file):
    """
    Create a plot showing sequence conservation across alignment positions.
    
    Args:
        alignment_file (str): Path to alignment file
        output_file (str): Path to save visualization
    """
    print("\nCreating conservation plot...")
    
    # Read alignment
    alignment = AlignIO.read(alignment_file, "fasta")
    
    alignment_length = alignment.get_alignment_length()
    num_sequences = len(alignment)
    
    # Calculate conservation score for each position
    conservation_scores = []
    
    for i in range(alignment_length):
        column = alignment[:, i]
        
        # Count non-gap residues
        residues = [aa for aa in column if aa != '-']
        
        if len(residues) == 0:
            conservation_scores.append(0)
            continue
        
        # Calculate conservation as frequency of most common residue
        from collections import Counter
        counts = Counter(residues)
        most_common_count = counts.most_common(1)[0][1]
        conservation = most_common_count / len(residues)
        conservation_scores.append(conservation)
    
    # Create plot
    fig, ax = plt.subplots(figsize=(20, 6))
    
    positions = range(1, alignment_length + 1)
    ax.plot(positions, conservation_scores, linewidth=0.5, color='blue')
    ax.fill_between(positions, conservation_scores, alpha=0.3, color='blue')
    
    ax.set_xlabel('Alignment Position', fontsize=12)
    ax.set_ylabel('Conservation Score', fontsize=12)
    ax.set_title('Sequence Conservation Across Alignment', fontsize=14, fontweight='bold')
    ax.set_ylim([0, 1])
    ax.grid(True, alpha=0.3)
    
    # Add horizontal line at 0.5
    ax.axhline(y=0.5, color='red', linestyle='--', alpha=0.5, label='50% conservation')
    ax.legend()
    
    plt.tight_layout()
    
    # Save figure
    conservation_plot_file = output_file.replace('.png', '_conservation.png')
    plt.savefig(conservation_plot_file, dpi=300, bbox_inches='tight')
    print(f"Conservation plot saved to: {conservation_plot_file}")
    
    plt.close()


def create_sequence_identity_matrix(alignment_file, output_file):
    """
    Create a heatmap of pairwise sequence identity.
    
    Args:
        alignment_file (str): Path to alignment file
        output_file (str): Path to save visualization
    """
    print("\nCreating sequence identity matrix...")
    
    # Read alignment
    alignment = AlignIO.read(alignment_file, "fasta")
    
    num_sequences = min(len(alignment), 50)  # Limit for readability
    
    # Calculate pairwise identity
    identity_matrix = np.zeros((num_sequences, num_sequences))
    
    for i in range(num_sequences):
        for j in range(num_sequences):
            seq_i = str(alignment[i].seq)
            seq_j = str(alignment[j].seq)
            
            matches = sum(1 for a, b in zip(seq_i, seq_j) if a == b and a != '-')
            total = sum(1 for a, b in zip(seq_i, seq_j) if a != '-' or b != '-')
            
            if total > 0:
                identity_matrix[i, j] = matches / total
            else:
                identity_matrix[i, j] = 0
    
    # Create plot
    fig, ax = plt.subplots(figsize=(12, 10))
    
    im = ax.imshow(identity_matrix, cmap='YlOrRd', vmin=0, vmax=1)
    
    ax.set_xlabel('Sequence Index', fontsize=12)
    ax.set_ylabel('Sequence Index', fontsize=12)
    ax.set_title('Pairwise Sequence Identity Matrix', fontsize=14, fontweight='bold')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Sequence Identity', fontsize=12)
    
    plt.tight_layout()
    
    # Save figure
    identity_plot_file = output_file.replace('.png', '_identity.png')
    plt.savefig(identity_plot_file, dpi=300, bbox_inches='tight')
    print(f"Identity matrix saved to: {identity_plot_file}")
    
    plt.close()


def main():
    """Main execution function."""
    # Check if alignment file exists
    if not os.path.exists(config.MSA_FILE):
        print(f"Error: MSA file not found: {config.MSA_FILE}")
        print("Please run generate_msa.py first")
        return
    
    # Define output file
    output_file = f"{config.RESULTS_DIR}/msa_heatmap.png"
    
    # Create visualizations
    visualize_msa_heatmap(config.MSA_FILE, output_file)
    create_conservation_plot(config.MSA_FILE, output_file)
    create_sequence_identity_matrix(config.MSA_FILE, output_file)
    
    print("\n" + "="*50)
    print("MSA visualization completed!")
    print(f"Visualizations saved in: {config.RESULTS_DIR}/")
    print("="*50)


if __name__ == "__main__":
    main()
