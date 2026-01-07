# Dehalogenases Investigation

A bioinformatics pipeline for investigating dehalogenase enzymes through sequence analysis and multiple sequence alignment.

## Overview

This repository provides tools to:
1. Fetch the protein sequence from PDB structure **3BJX** (a dehalogenase enzyme)
2. Perform BLAST searches to identify homologous sequences from protein databases
3. Generate multiple sequence alignments (MSA) using Clustal Omega
4. Visualize the MSA with various plots and heatmaps

## Requirements

### Python Dependencies
Install the required Python packages:

```bash
pip install -r requirements.txt
```

Required packages:
- **biopython**: For sequence manipulation and BLAST
- **matplotlib**: For visualization
- **numpy**: For numerical operations
- **requests**: For fetching PDB data

### External Tools
For proper MSA generation, install Clustal Omega:

**Ubuntu/Debian:**
```bash
sudo apt-get install clustalo
```

**macOS:**
```bash
brew install clustal-omega
```

**Note:** The pipeline will work without Clustal Omega but will produce unaligned sequences.

## Project Structure

```
dehalogenases/
├── src/
│   ├── fetch_pdb_sequence.py    # Fetch sequence from PDB
│   ├── blast_search.py           # Run BLAST searches
│   ├── generate_msa.py           # Generate MSA
│   └── visualize_msa.py          # Create visualizations
├── data/                         # Raw data (sequences, BLAST results)
├── results/                      # Output files (alignments, plots)
├── config.py                     # Configuration parameters
├── run_pipeline.py              # Main pipeline script
├── requirements.txt             # Python dependencies
└── README.md                    # This file
```

## Usage

### Quick Start: Run Complete Pipeline

Run the entire pipeline with one command:

```bash
python run_pipeline.py
```

**Note:** The BLAST search step can take 5-10 minutes. To skip it during testing:

```bash
python run_pipeline.py --skip-blast
```

### Step-by-Step Execution

You can also run each step individually:

1. **Fetch PDB Sequence:**
   ```bash
   python src/fetch_pdb_sequence.py
   ```

2. **Run BLAST Search:**
   ```bash
   python src/blast_search.py
   ```
   *Note: This step takes several minutes*

3. **Generate MSA:**
   ```bash
   python src/generate_msa.py
   ```

4. **Visualize MSA:**
   ```bash
   python src/visualize_msa.py
   ```

## Output Files

After running the pipeline, you'll find:

### Data Directory (`data/`)
- `3BJX_sequence.fasta`: Original PDB sequence
- `blast_results.xml`: Raw BLAST search results
- `homologs.fasta`: Homologous sequences from BLAST
- `combined_sequences.fasta`: All sequences for alignment

### Results Directory (`results/`)
- `alignment.fasta`: Multiple sequence alignment (MSA)
- `msa_heatmap.png`: Heatmap visualization of the alignment
- `msa_heatmap_conservation.png`: Conservation plot showing conserved regions
- `msa_heatmap_identity.png`: Pairwise sequence identity matrix

## Configuration

Edit `config.py` to customize:
- PDB structure ID
- BLAST parameters (e-value, max hits, database)
- MSA tool settings
- Output file paths

## About Dehalogenases

Dehalogenases are enzymes that catalyze the removal of halogen atoms from organic compounds. They are important for:
- Biodegradation of environmental pollutants
- Industrial biocatalysis
- Understanding enzyme evolution and mechanism

**PDB 3BJX** is a haloalkane dehalogenase structure that serves as the reference sequence for this analysis.

## Visualization Features

The pipeline generates three types of visualizations:

1. **MSA Heatmap**: Color-coded alignment showing:
   - Hydrophobic residues (orange)
   - Polar residues (green)
   - Positively charged residues (blue)
   - Negatively charged residues (red)
   - Gaps (white)

2. **Conservation Plot**: Shows sequence conservation across alignment positions

3. **Identity Matrix**: Heatmap of pairwise sequence identity between all sequences

## Troubleshooting

### BLAST Search Takes Too Long
- The BLAST search against UniProt can take 5-10 minutes
- Use `--skip-blast` flag to test other parts of the pipeline
- Reduce `MAX_HITS` in `config.py` for faster searches

### Clustal Omega Not Found
- Install Clustal Omega (see Requirements section)
- The pipeline will create unaligned sequences if not installed

### Memory Issues
- Large alignments may require substantial memory
- Reduce `MAX_HITS` in `config.py`
- Visualizations automatically limit display size

## License

This project is provided for educational and research purposes.

## Author

Created for investigating dehalogenase enzyme evolution and sequence diversity.
