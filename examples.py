"""
Example usage of the dehalogenases package.

This script demonstrates basic functionality of the package for analyzing
dehalogenase enzyme sequences.
"""

from dehalogenases import Sequence, SequenceAnalyzer


def main():
    """Demonstrate basic usage of the dehalogenases package."""

    print("=" * 60)
    print("Dehalogenases Package Example")
    print("=" * 60)
    print()

    # Create some example dehalogenase sequences
    sequences = [
        Sequence(
            id="HAD_ECOLI_001",
            sequence="MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIHHLAHSLVGKRRRPLTIDSSIEPDYKKDDDDKAM",
            organism="Escherichia coli",
            enzyme_class="HAD",
        ),
        Sequence(
            id="HAD_BSUBT_001",
            sequence="MKRLIEASSLREFKPQSAFELIAANADGSAIDLSYQERIKEGVAIDGTTTATVLQAVKMGAEVDIVVFDDADILDEFERQQGRTNLTSAIHKLGIGIHRSGESRTTYYRSLLEDIQPGDGKKFKNMKLSFESAAGGEYDVPAEDAARVYNKGLKYVVAGDAPAGKNKKGLLHIAGKDLRSRLKERL",
            organism="Bacillus subtilis",
            enzyme_class="HAD",
        ),
        Sequence(
            id="SC_PAERU_001",
            sequence="MTQPKPSNLAAARLFGLAAVAGTALSGCSGKEVYDYDASQGEKVAIVGSDADLAEALANGASVIVSGGSAMTTDAIADAARNGGTKVLMPGSTSREGETAMVVDKGTAKAVVMDFSNADWQELRPALKQ",
            organism="Pseudomonas aeruginosa",
            enzyme_class="short-chain",
        ),
        Sequence(
            id="HAD_PAERU_002",
            sequence="MKLGIAVGLAAAVMPQLAGCAPQQAYEVFDASRDFVMVVGDDGFIDKEREAAGLHVEMGGSRTRTTAAEVDAGGTRLVLPDSATYEGEEALTYDASVATATVFEGDKLHQSLREMKDALQKP",
            organism="Pseudomonas aeruginosa",
            enzyme_class="HAD",
        ),
    ]

    # Demonstrate individual sequence analysis
    print("Individual Sequence Analysis:")
    print("-" * 60)
    seq = sequences[0]
    print(f"ID: {seq.id}")
    print(f"Organism: {seq.organism}")
    print(f"Class: {seq.enzyme_class}")
    print(f"Length: {len(seq)} amino acids")
    print(f"First 40 residues: {seq.sequence[:40]}...")
    print()

    composition = seq.get_composition()
    print(f"Most common amino acids:")
    sorted_aa = sorted(composition.items(), key=lambda x: x[1], reverse=True)[:5]
    for aa, count in sorted_aa:
        print(f"  {aa}: {count} ({count/len(seq)*100:.1f}%)")
    print()

    # Demonstrate sequence analyzer
    print("Multi-Sequence Analysis:")
    print("-" * 60)
    analyzer = SequenceAnalyzer(sequences=sequences)

    stats = analyzer.get_length_statistics()
    print(f"Total sequences: {stats['count']}")
    print(f"Length range: {stats['min']} - {stats['max']} amino acids")
    print(f"Mean length: {stats['mean']:.1f} amino acids")
    print()

    organisms = analyzer.get_unique_organisms()
    print(f"Organisms represented ({len(organisms)}):")
    for org in organisms:
        count = len(analyzer.filter_by_organism(org))
        print(f"  - {org}: {count} sequence(s)")
    print()

    # Filter by enzyme class
    print("Filtering by Enzyme Class:")
    print("-" * 60)
    had_sequences = analyzer.filter_by_class("HAD")
    print(f"HAD family sequences: {len(had_sequences)}")
    for seq in had_sequences:
        print(f"  - {seq.id}: {seq.organism} ({len(seq)} aa)")
    print()

    sc_sequences = analyzer.filter_by_class("short-chain")
    print(f"Short-chain family sequences: {len(sc_sequences)}")
    for seq in sc_sequences:
        print(f"  - {seq.id}: {seq.organism} ({len(seq)} aa)")
    print()

    print("=" * 60)
    print("Example completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()
