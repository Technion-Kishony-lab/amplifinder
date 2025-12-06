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


def make_tn_sequences() -> dict:
    """Define TN element sequences for testing."""
    return {
        "IS_test1": make_random_seq(700, seed=100),
        "IS_test2": make_random_seq(600, seed=101),
    }


def make_reference_with_tn(tn_seqs: dict) -> tuple[Seq, dict]:
    """Build reference genome with embedded TN elements.

    Returns:
        (sequence, tn_locations) where tn_locations maps TN name to (start, end, complement)
    """
    # Flanking regions
    flank1 = make_random_seq(500, seed=1)
    flank2 = make_random_seq(400, seed=2)
    flank3 = make_random_seq(500, seed=3)

    tn1 = tn_seqs["IS_test1"]  # forward strand
    tn2_rc = tn_seqs["IS_test2"].reverse_complement()

    # Build sequence and track positions (1-based)
    seq = flank1 + tn1 + flank2 + tn2_rc + flank3

    tn_locations = {
        "IS_test1": (501, 500 + len(tn1), False),           # forward
        "IS_test2": (501 + len(tn1) + 400, 500 + len(tn1) + 400 + len(tn2_rc), True),  # complement
    }

    return seq, tn_locations


def write_fasta(path: Path, name: str, seq: Seq) -> None:
    """Write single-sequence FASTA file using BioPython."""
    record = SeqRecord(seq, id=name, description="")
    SeqIO.write(record, path, "fasta")


def write_TN_database(path: Path, tn_seqs: dict) -> None:
    """Write TN database FASTA using BioPython."""
    records = [SeqRecord(seq, id=name, description="") for name, seq in tn_seqs.items()]
    SeqIO.write(records, path, "fasta")


def write_genbank(path: Path, name: str, seq: Seq, tn_locations: dict) -> None:
    """Write GenBank file with TN annotations using BioPython."""
    record = SeqRecord(
        seq,
        id=name,
        name=name,
        description="Synthetic reference with TN elements for testing",
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

    # Add TN features
    for tn_name, (start, end, complement) in tn_locations.items():
        strand = -1 if complement else 1
        record.features.append(SeqFeature(
            FeatureLocation(start - 1, end, strand=strand),  # BioPython uses 0-based
            type="mobile_element",
            qualifiers={
                "mobile_element_type": [f"insertion sequence:{tn_name}"],
                "note": ["Test TN element"],
            },
        ))

    SeqIO.write(record, path, "genbank")


# =============================================================================
# Fixtures
# =============================================================================

@pytest.fixture
def tn_sequences():
    """TN element sequences used in tests."""
    return make_tn_sequences()


@pytest.fixture
def test_data(tmp_path, tn_sequences):
    """Generate synthetic test data in tmp_path.

    Creates:
        - tiny_ref.fasta: Reference with embedded TN elements
        - tiny_ref.gbk: GenBank with TN annotations
        - tiny_tn.fna: TN database

    Returns:
        Path to test_data directory
    """
    data_dir = tmp_path / "test_data"
    data_dir.mkdir()

    # Build reference with TN elements
    seq, tn_locs = make_reference_with_tn(tn_sequences)

    # Write files (use "tiny" as sequence ID to match ref_name in tests)
    write_fasta(data_dir / "tiny_ref.fasta", "tiny", seq)
    write_genbank(data_dir / "tiny_ref.gbk", "tiny", seq, tn_locs)
    write_TN_database(data_dir / "tiny_tn.fna", tn_sequences)

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
def tiny_tn_db(test_data):
    """Path to synthetic TN database."""
    return test_data / "tiny_tn.fna"
