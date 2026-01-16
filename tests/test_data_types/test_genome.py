import pytest
from Bio.Seq import reverse_complement

from amplifinder.data_types.genome import Genome
from amplifinder.data_types.enums import Orientation
from amplifinder.data_types.record_types import Junction


@pytest.fixture
def toy_genome(tmp_path) -> Genome:
    """Create a toy genome for testing."""
    seq = "".join(chr(ord('a') + i) for i in range(26))  # abc...xyz
    fasta = tmp_path / "toy.fasta"
    fasta.write_text(f">toy\n{seq}\n")
    return Genome(name="toy", fasta_path=fasta)


def test_get_forward_sequence_in_range(toy_genome: Genome):
    scaf = toy_genome.get_scaffold("toy")
    assert scaf.slice(3, 8) == "cdefgh"


def test_get_forward_sequence_in_range_circular(toy_genome: Genome):
    scaf = toy_genome.get_scaffold("toy")
    assert scaf.slice(25, 2) == "yzab"


def test_get_sequence_in_range_reverse(toy_genome: Genome):
    scaf = toy_genome.get_scaffold("toy")
    result = scaf.slice(10, 5, Orientation.REVERSE)
    assert result == reverse_complement("efghij")


def test_junction_sequences(toy_genome: Genome):
    jc = Junction(
        scaf1="toy",
        pos1=10,
        dir1=Orientation.REVERSE,
        scaf2="toy",
        pos2=11,
        dir2=Orientation.FORWARD,
        flanking1=5,
        flanking2=5,
    )

    arm1 = toy_genome.get_junction_arm_sequence(jc, 1)
    arm2 = toy_genome.get_junction_arm_sequence(jc, 2)

    assert arm1 == reverse_complement("fghij")
    assert arm2 == "klmno"
    assert toy_genome.get_junction_sequence_arm1_to_arm2(jc) == "fghijklmno"
    assert toy_genome.get_junction_sequence_arm2_to_arm1(jc) == reverse_complement("fghijklmno")
