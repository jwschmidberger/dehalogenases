"""
Dehalogenases: A collection of code for evolutionary investigations into dehalogenase enzymes.

This package provides tools for analyzing, comparing, and studying the evolution
of dehalogenase enzymes across different organisms.
"""

__version__ = "0.1.0"

from .sequence import Sequence
from .analysis import SequenceAnalyzer

__all__ = ["Sequence", "SequenceAnalyzer", "__version__"]
