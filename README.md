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
