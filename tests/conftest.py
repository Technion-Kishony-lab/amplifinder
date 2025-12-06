"""Pytest fixtures for AmpliFinder tests."""

import random
import pytest
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

# =============================================================================
# Assertion helpers
# =============================================================================


# =============================================================================
# Test configuration
# =============================================================================

FIXTURES_DIR = Path(__file__).parent / "fixtures"


# =============================================================================
# Synthetic test data generators
# =============================================================================

def make_random_seq(length: int, seed: int = 42) -> Seq:
    """Generate reproducible random DNA sequence."""
    rng = random.Random(seed)
    return Seq("".join(rng.choices("ACGT", k=length)))


def make_is_sequences() -> dict:
    """Define IS element sequences for testing."""
    return {
        "IS_test1": make_random_seq(700, seed=100),
        "IS_test2": make_random_seq(600, seed=101),
    }


def make_reference_with_is(is_seqs: dict) -> tuple[Seq, dict]:
    """Build reference genome with embedded IS elements.

    Returns:
        (sequence, is_locations) where is_locations maps IS name to (start, end, complement)
    """
    # Flanking regions
    flank1 = make_random_seq(500, seed=1)
    flank2 = make_random_seq(400, seed=2)
    flank3 = make_random_seq(500, seed=3)

    is1 = is_seqs["IS_test1"]  # forward strand
    is2_rc = is_seqs["IS_test2"].reverse_complement()

    # Build sequence and track positions (1-based)
    seq = flank1 + is1 + flank2 + is2_rc + flank3

    is_locations = {
        "IS_test1": (501, 500 + len(is1), False),           # forward
        "IS_test2": (501 + len(is1) + 400, 500 + len(is1) + 400 + len(is2_rc), True),  # complement
    }

    return seq, is_locations


def write_fasta(path: Path, name: str, seq: Seq) -> None:
    """Write single-sequence FASTA file using BioPython."""
    record = SeqRecord(seq, id=name, description="")
    SeqIO.write(record, path, "fasta")


def write_IS_database(path: Path, is_seqs: dict) -> None:
    """Write IS database FASTA using BioPython."""
    records = [SeqRecord(seq, id=name, description="") for name, seq in is_seqs.items()]
    SeqIO.write(records, path, "fasta")


def write_genbank(path: Path, name: str, seq: Seq, is_locations: dict) -> None:
    """Write GenBank file with IS annotations using BioPython."""
    record = SeqRecord(
        seq,
        id=name,
        name=name,
        description="Synthetic reference with IS elements for testing",
        annotations={
            "molecule_type": "DNA",
            "topology": "linear",
            "data_file_division": "SYN",
        },
    )

    # Add source feature
    record.features.append(SeqFeature(
        FeatureLocation(0, len(seq)),
        type="source",
        qualifiers={"organism": ["synthetic construct"], "mol_type": ["genomic DNA"]},
    ))

    # Add IS features
    for is_name, (start, end, complement) in is_locations.items():
        strand = -1 if complement else 1
        record.features.append(SeqFeature(
            FeatureLocation(start - 1, end, strand=strand),  # BioPython uses 0-based
            type="mobile_element",
            qualifiers={
                "mobile_element_type": [f"insertion sequence:{is_name}"],
                "note": ["Test IS element"],
            },
        ))

    SeqIO.write(record, path, "genbank")


# =============================================================================
# Fixtures
# =============================================================================

@pytest.fixture
def is_sequences():
    """IS element sequences used in tests."""
    return make_is_sequences()


@pytest.fixture
def test_data(tmp_path, is_sequences):
    """Generate synthetic test data in tmp_path.

    Creates:
        - tiny_ref.fasta: Reference with embedded IS elements
        - tiny_ref.gbk: GenBank with IS annotations
        - tiny_is.fna: IS database

    Returns:
        Path to test_data directory
    """
    data_dir = tmp_path / "test_data"
    data_dir.mkdir()

    # Build reference with IS elements
    seq, is_locs = make_reference_with_is(is_sequences)

    # Write files (use "tiny" as sequence ID to match ref_name in tests)
    write_fasta(data_dir / "tiny_ref.fasta", "tiny", seq)
    write_genbank(data_dir / "tiny_ref.gbk", "tiny", seq, is_locs)
    write_IS_database(data_dir / "tiny_is.fna", is_sequences)

    return data_dir


@pytest.fixture
def tmp_output(tmp_path):
    """Isolated output directory for each test."""
    out = tmp_path / "output"
    out.mkdir()
    return out


@pytest.fixture
def tiny_ref_fasta(test_data):
    """Path to synthetic reference FASTA."""
    return test_data / "tiny_ref.fasta"


@pytest.fixture
def tiny_ref_gbk(test_data):
    """Path to synthetic reference GenBank."""
    return test_data / "tiny_ref.gbk"


@pytest.fixture
def tiny_is_db(test_data):
    """Path to synthetic IS database."""
    return test_data / "tiny_is.fna"
