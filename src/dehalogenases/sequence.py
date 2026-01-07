"""
Sequence module for handling dehalogenase enzyme sequences.
"""

from typing import Optional


class Sequence:
    """
    Represents a dehalogenase enzyme sequence.

    Attributes:
        id: Unique identifier for the sequence
        sequence: The amino acid sequence string
        organism: Optional organism name
        enzyme_class: Optional enzyme classification
    """

    def __init__(
        self,
        id: str,
        sequence: str,
        organism: Optional[str] = None,
        enzyme_class: Optional[str] = None,
    ):
        """
        Initialize a Sequence object.

        Args:
            id: Unique identifier for the sequence
            sequence: The amino acid sequence string
            organism: Optional organism name
            enzyme_class: Optional enzyme classification (e.g., 'HAD', 'short-chain')
        """
        self.id = id
        self.sequence = sequence.upper()
        self.organism = organism
        self.enzyme_class = enzyme_class

    def __len__(self) -> int:
        """Return the length of the sequence."""
        return len(self.sequence)

    def __str__(self) -> str:
        """Return a string representation of the sequence."""
        return f"Sequence(id={self.id}, length={len(self)}, organism={self.organism})"

    def __repr__(self) -> str:
        """Return a detailed representation of the sequence."""
        return (
            f"Sequence(id='{self.id}', sequence='{self.sequence[:20]}...', "
            f"organism='{self.organism}', enzyme_class='{self.enzyme_class}')"
        )

    def get_composition(self) -> dict[str, int]:
        """
        Calculate amino acid composition of the sequence.

        Returns:
            Dictionary mapping amino acid codes to their counts
        """
        composition = {}
        for aa in self.sequence:
            composition[aa] = composition.get(aa, 0) + 1
        return composition

    def get_gc_residue_proportion(self) -> float:
        """
        Calculate the proportion of glycine (G) and cysteine (C) residues.

        This can be useful for analyzing protein properties, as glycine provides
        flexibility and cysteine is important for disulfide bond formation.

        Returns:
            Proportion of G and C residues (0.0 to 1.0)
        """
        gc_count = self.sequence.count("G") + self.sequence.count("C")
        if len(self.sequence) == 0:
            return 0.0
        return gc_count / len(self.sequence)
