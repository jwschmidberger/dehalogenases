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
# Dehalogenases

A collection of code focused on evolutionary investigations into dehalogenase enzymes.

## Overview

This package provides tools for analyzing, comparing, and studying the evolution of dehalogenase enzymes across different organisms. Dehalogenases are enzymes that catalyze the removal of halogen atoms from organic compounds, playing important roles in bioremediation and industrial processes.

## Installation

### Development Installation

Clone the repository and install in development mode:

```bash
git clone https://github.com/jwschmidberger/dehalogenases.git
cd dehalogenases
pip install -e ".[dev]"
```

### Basic Installation

```bash
pip install -e .
```

## Usage

### Basic Example

```python
from dehalogenases import Sequence, SequenceAnalyzer

# Create a dehalogenase sequence
seq1 = Sequence(
    id="DEH001",
    sequence="MKTAYIAKQRQISFVKSHFSRQ",
    organism="E. coli",
    enzyme_class="HAD"
)

# Analyze the sequence
print(f"Length: {len(seq1)}")
print(f"Composition: {seq1.get_composition()}")

# Create an analyzer with multiple sequences
seq2 = Sequence(
    id="DEH002",
    sequence="YIAKQRISFVKSHFSRQMKTA",
    organism="B. subtilis",
    enzyme_class="HAD"
)

analyzer = SequenceAnalyzer(sequences=[seq1, seq2])
print(f"Average length: {analyzer.get_average_length()}")
print(f"Statistics: {analyzer.get_length_statistics()}")
print(f"Organisms: {analyzer.get_unique_organisms()}")
```

## Development

### Running Tests

```bash
pytest
```

### Code Formatting

```bash
black src tests
```

### Linting

```bash
ruff check src tests
```

## Project Structure

```
dehalogenases/
├── src/
│   └── dehalogenases/
│       ├── __init__.py
│       ├── sequence.py      # Sequence representation and basic operations
│       └── analysis.py      # Analysis tools for sequence collections
├── tests/
│   ├── test_sequence.py
│   └── test_analysis.py
├── pyproject.toml           # Project configuration and dependencies
└── README.md
```

## Features

- **Sequence Management**: Store and manipulate dehalogenase enzyme sequences
- **Composition Analysis**: Calculate amino acid composition and properties
- **Comparative Analysis**: Analyze multiple sequences together
- **Filtering**: Filter sequences by organism or enzyme class
- **Statistics**: Calculate length statistics and distributions

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is open source and available under the MIT License.
