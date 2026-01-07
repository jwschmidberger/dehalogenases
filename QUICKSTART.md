# Quick Start Guide

This guide will help you get started with the dehalogenases investigation pipeline.

## Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/jwschmidberger/dehalogenases.git
   cd dehalogenases
   ```

2. **Install Python dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

3. **(Optional) Install Clustal Omega for proper sequence alignment:**
   ```bash
   # Ubuntu/Debian
   sudo apt-get install clustalo
   
   # macOS
   brew install clustal-omega
   ```

## Quick Demo

Run the example script to see the complete workflow:

```bash
python example.py
```

This will:
1. Fetch the 3BJX sequence (haloalkane dehalogenase)
2. Use demo homologous sequences
3. Generate a multiple sequence alignment
4. Create visualizations

## Full Pipeline

To run the complete pipeline with all steps:

```bash
python run_pipeline.py
```

**Note:** The BLAST search step can take 5-10 minutes. To skip it:
```bash
python run_pipeline.py --skip-blast
```

## Individual Scripts

Run each step separately:

### 1. Fetch PDB Sequence
```bash
python src/fetch_pdb_sequence.py
```

Fetches the protein sequence for PDB structure 3BJX and saves it to `data/3BJX_sequence.fasta`.

### 2. BLAST Search (Optional)
```bash
python src/blast_search.py
```

Searches for homologous sequences in UniProt. This step requires internet access and takes 5-10 minutes.

**For demo purposes, homologous sequences are already provided in `data/homologs.fasta`.**

### 3. Generate Multiple Sequence Alignment
```bash
python src/generate_msa.py
```

Creates a multiple sequence alignment from the query and homologous sequences.

### 4. Visualize MSA
```bash
python src/visualize_msa.py
```

Generates three visualization files:
- `results/msa_heatmap.png` - Color-coded alignment heatmap
- `results/msa_heatmap_conservation.png` - Conservation plot
- `results/msa_heatmap_identity.png` - Pairwise identity matrix

## Output Files

### Data Directory (`data/`)
- `3BJX_sequence.fasta` - Query sequence from PDB
- `homologs.fasta` - Homologous dehalogenase sequences
- `combined_sequences.fasta` - All sequences for alignment
- `blast_results.xml` - Raw BLAST results (if BLAST search was run)

### Results Directory (`results/`)
- `alignment.fasta` - Multiple sequence alignment
- `msa_heatmap.png` - Alignment visualization
- `msa_heatmap_conservation.png` - Conservation analysis
- `msa_heatmap_identity.png` - Sequence similarity matrix

## Customization

Edit `config.py` to customize:

```python
# PDB structure to analyze
PDB_ID = "3BJX"

# BLAST parameters
BLAST_DATABASE = "swissprot"
E_VALUE_THRESHOLD = 0.001
MAX_HITS = 100

# Output directories
DATA_DIR = "data"
RESULTS_DIR = "results"
```

## Interpreting Results

### MSA Heatmap
- **Orange**: Hydrophobic residues (A, V, L, I, M, F, W, P)
- **Green**: Polar residues (S, T, N, Q, C, G, Y)
- **Blue**: Positively charged (K, R, H)
- **Red**: Negatively charged (D, E)
- **White**: Gaps in alignment
- **Gray**: Unknown/ambiguous residues

### Conservation Plot
- Shows conservation score (0-1) for each alignment position
- Higher values = more conserved across sequences
- Red dashed line indicates 50% conservation threshold
- Highly conserved regions often represent functionally important sites

### Identity Matrix
- Heatmap showing pairwise sequence identity
- Darker red = higher identity between sequences
- Diagonal is always 1.0 (100% identity with itself)
- Helps identify closely related sequences

## Troubleshooting

### "Clustal Omega not found"
The pipeline will work without Clustal Omega but will create a simple demo alignment. For production use, install Clustal Omega as shown in the Installation section.

### "No internet connection"
The `fetch_pdb_sequence.py` script automatically falls back to an offline mode with pre-loaded 3BJX sequence. BLAST search requires internet access but demo data is provided.

### Memory issues with large alignments
Reduce `MAX_HITS` in `config.py` to limit the number of sequences in the alignment.

## Next Steps

1. Explore the alignment file to identify conserved regions
2. Analyze the visualizations to understand sequence diversity
3. Run BLAST search with your own parameters
4. Try different PDB structures (modify `config.py`)
5. Extend the analysis with phylogenetic tree construction

## Support

For issues or questions, please open an issue on GitHub.
