import pytest
from Bio.Seq import reverse_complement

from amplifinder.data_types.genome import Genome
from amplifinder.data_types.record_types import Junction, Orientation


@pytest.fixture
def toy_genome(tmp_path) -> Genome:
    """Create a toy genome for testing."""
    seq = "".join(chr(ord('a') + i) for i in range(26))  # abc...xyz
    fasta = tmp_path / "toy.fasta"
    fasta.write_text(f">toy\n{seq}\n")
    return Genome(name="toy", fasta_path=fasta)


def test_get_forward_sequence_in_range(toy_genome):
    assert toy_genome.get_fowrard_sequence_in_range("toy", 3, 8) == "cdefgh"


def test_get_forward_sequence_in_range_circular(toy_genome):
    assert toy_genome.get_fowrard_sequence_in_range("toy", 25, 2) == "yzab"


def test_get_sequence_in_range_reverse(toy_genome):
    result = toy_genome.get_sequence_in_range("toy", 10, 5, Orientation.REVERSE)
    print(reverse_complement("efghij"))
    assert result == reverse_complement("efghij")


def test_junction_sequences(toy_genome):
    jc = Junction(
        num=0,
        scaf1="toy",
        pos1=10,
        dir1=Orientation.REVERSE,
        scaf2="toy",
        pos2=11,
        dir2=Orientation.FORWARD,
        flanking_left=5,
        flanking_right=5,
    )

    side1 = toy_genome.get_junction_side_sequence(jc, 1)
    side2 = toy_genome.get_junction_side_sequence(jc, 2)

    assert side1 == reverse_complement("efghij")
    assert side2 == "klmnop"
    assert toy_genome.get_junction_sequence(jc) == "efghijklmnop"

