"""
Tests for the analysis module.
"""

from dehalogenases.sequence import Sequence
from dehalogenases.analysis import SequenceAnalyzer


class TestSequenceAnalyzer:
    """Test cases for the SequenceAnalyzer class."""

    def test_analyzer_initialization_empty(self):
        """Test analyzer initialization with no sequences."""
        analyzer = SequenceAnalyzer()
        assert analyzer.sequences == []

    def test_analyzer_initialization_with_sequences(self):
        """Test analyzer initialization with sequences."""
        seq1 = Sequence(id="DEH001", sequence="MKTA")
        seq2 = Sequence(id="DEH002", sequence="YIAK")
        analyzer = SequenceAnalyzer(sequences=[seq1, seq2])
        assert len(analyzer.sequences) == 2

    def test_add_sequence(self):
        """Test adding a sequence to analyzer."""
        analyzer = SequenceAnalyzer()
        seq = Sequence(id="DEH003", sequence="MKTA")
        analyzer.add_sequence(seq)
        assert len(analyzer.sequences) == 1
        assert analyzer.sequences[0].id == "DEH003"

    def test_get_average_length(self):
        """Test average length calculation."""
        seq1 = Sequence(id="DEH004", sequence="MKTA")  # length 4
        seq2 = Sequence(id="DEH005", sequence="YIAKQR")  # length 6
        analyzer = SequenceAnalyzer(sequences=[seq1, seq2])
        assert analyzer.get_average_length() == 5.0

    def test_get_average_length_empty(self):
        """Test average length with no sequences."""
        analyzer = SequenceAnalyzer()
        assert analyzer.get_average_length() == 0.0

    def test_get_length_statistics(self):
        """Test length statistics calculation."""
        seq1 = Sequence(id="DEH006", sequence="MKTA")  # length 4
        seq2 = Sequence(id="DEH007", sequence="YIAKQR")  # length 6
        seq3 = Sequence(id="DEH008", sequence="YI")  # length 2
        analyzer = SequenceAnalyzer(sequences=[seq1, seq2, seq3])

        stats = analyzer.get_length_statistics()
        assert stats["count"] == 3
        assert stats["min"] == 2
        assert stats["max"] == 6
        assert stats["mean"] == 4.0

    def test_get_length_statistics_empty(self):
        """Test length statistics with no sequences."""
        analyzer = SequenceAnalyzer()
        stats = analyzer.get_length_statistics()
        assert stats["count"] == 0
        assert stats["min"] == 0
        assert stats["max"] == 0
        assert stats["mean"] == 0.0

    def test_filter_by_organism(self):
        """Test filtering sequences by organism."""
        seq1 = Sequence(id="DEH009", sequence="MKTA", organism="E. coli")
        seq2 = Sequence(id="DEH010", sequence="YIAK", organism="B. subtilis")
        seq3 = Sequence(id="DEH011", sequence="QRSF", organism="E. coli")
        analyzer = SequenceAnalyzer(sequences=[seq1, seq2, seq3])

        ecoli_seqs = analyzer.filter_by_organism("E. coli")
        assert len(ecoli_seqs) == 2
        assert all(seq.organism == "E. coli" for seq in ecoli_seqs)

    def test_filter_by_class(self):
        """Test filtering sequences by enzyme class."""
        seq1 = Sequence(id="DEH012", sequence="MKTA", enzyme_class="HAD")
        seq2 = Sequence(id="DEH013", sequence="YIAK", enzyme_class="short-chain")
        seq3 = Sequence(id="DEH014", sequence="QRSF", enzyme_class="HAD")
        analyzer = SequenceAnalyzer(sequences=[seq1, seq2, seq3])

        had_seqs = analyzer.filter_by_class("HAD")
        assert len(had_seqs) == 2
        assert all(seq.enzyme_class == "HAD" for seq in had_seqs)

    def test_get_unique_organisms(self):
        """Test getting unique organism names."""
        seq1 = Sequence(id="DEH015", sequence="MKTA", organism="E. coli")
        seq2 = Sequence(id="DEH016", sequence="YIAK", organism="B. subtilis")
        seq3 = Sequence(id="DEH017", sequence="QRSF", organism="E. coli")
        seq4 = Sequence(id="DEH018", sequence="VKSH", organism="P. aeruginosa")
        analyzer = SequenceAnalyzer(sequences=[seq1, seq2, seq3, seq4])

        organisms = analyzer.get_unique_organisms()
        assert len(organisms) == 3
        assert "E. coli" in organisms
        assert "B. subtilis" in organisms
        assert "P. aeruginosa" in organisms
        assert organisms == sorted(organisms)  # Check if sorted
