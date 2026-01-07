"""
Tests for the sequence module.
"""

from dehalogenases.sequence import Sequence


class TestSequence:
    """Test cases for the Sequence class."""

    def test_sequence_initialization(self):
        """Test basic sequence initialization."""
        seq = Sequence(
            id="DEH001",
            sequence="MKTAYIAKQRQISFVKSHFSRQ",
            organism="E. coli",
            enzyme_class="HAD",
        )
        assert seq.id == "DEH001"
        assert seq.sequence == "MKTAYIAKQRQISFVKSHFSRQ"
        assert seq.organism == "E. coli"
        assert seq.enzyme_class == "HAD"

    def test_sequence_uppercase_conversion(self):
        """Test that sequences are converted to uppercase."""
        seq = Sequence(id="DEH002", sequence="mktayiakqrq")
        assert seq.sequence == "MKTAYIAKQRQ"

    def test_sequence_length(self):
        """Test sequence length calculation."""
        seq = Sequence(id="DEH003", sequence="MKTAYIAKQRQISFVKSHFSRQ")
        assert len(seq) == 22

    def test_sequence_str(self):
        """Test string representation."""
        seq = Sequence(id="DEH004", sequence="MKTA", organism="E. coli")
        str_repr = str(seq)
        assert "DEH004" in str_repr
        assert "4" in str_repr  # length
        assert "E. coli" in str_repr

    def test_get_composition(self):
        """Test amino acid composition calculation."""
        seq = Sequence(id="DEH005", sequence="MKTAYK")
        composition = seq.get_composition()
        assert composition["M"] == 1
        assert composition["K"] == 2
        assert composition["T"] == 1
        assert composition["A"] == 1
        assert composition["Y"] == 1

    def test_gc_content(self):
        """Test GC content calculation."""
        seq = Sequence(id="DEH006", sequence="GGCC")
        assert seq.gc_content() == 1.0

        seq2 = Sequence(id="DEH007", sequence="GGAA")
        assert seq2.gc_content() == 0.5

        seq3 = Sequence(id="DEH008", sequence="AAAA")
        assert seq3.gc_content() == 0.0

    def test_empty_sequence(self):
        """Test handling of empty sequence."""
        seq = Sequence(id="DEH009", sequence="")
        assert len(seq) == 0
        assert seq.gc_content() == 0.0
        assert seq.get_composition() == {}
