"""
Analysis module for dehalogenase sequence analysis and comparison.
"""

from typing import Dict, List, Optional, Union
from .sequence import Sequence


class SequenceAnalyzer:
    """
    Provides analysis tools for dehalogenase sequences.
    """

    def __init__(self, sequences: Optional[List[Sequence]] = None):
        """
        Initialize a SequenceAnalyzer.

        Args:
            sequences: Optional list of Sequence objects to analyze
        """
        self.sequences = sequences if sequences is not None else []

    def add_sequence(self, sequence: Sequence) -> None:
        """
        Add a sequence to the analyzer.

        Args:
            sequence: Sequence object to add
        """
        self.sequences.append(sequence)

    def get_average_length(self) -> float:
        """
        Calculate the average length of sequences.

        Returns:
            Average sequence length
        """
        if not self.sequences:
            return 0.0
        return sum(len(seq) for seq in self.sequences) / len(self.sequences)

    def get_length_statistics(self) -> Dict[str, Union[int, float]]:
        """
        Get length statistics for all sequences.

        Returns:
            Dictionary containing min, max, mean, and count statistics
        """
        if not self.sequences:
            return {"count": 0, "min": 0, "max": 0, "mean": 0.0}

        lengths = [len(seq) for seq in self.sequences]
        return {
            "count": len(lengths),
            "min": min(lengths),
            "max": max(lengths),
            "mean": sum(lengths) / len(lengths),
        }

    def filter_by_organism(self, organism: str) -> List[Sequence]:
        """
        Filter sequences by organism name.

        Args:
            organism: Organism name to filter by

        Returns:
            List of sequences from the specified organism
        """
        return [seq for seq in self.sequences if seq.organism == organism]

    def filter_by_class(self, enzyme_class: str) -> List[Sequence]:
        """
        Filter sequences by enzyme class.

        Args:
            enzyme_class: Enzyme class to filter by

        Returns:
            List of sequences of the specified class
        """
        return [seq for seq in self.sequences if seq.enzyme_class == enzyme_class]

    def get_unique_organisms(self) -> List[str]:
        """
        Get a list of unique organism names in the dataset.

        Returns:
            List of unique organism names
        """
        organisms = set()
        for seq in self.sequences:
            if seq.organism:
                organisms.add(seq.organism)
        return sorted(list(organisms))
