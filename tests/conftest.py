"""Pytest fixtures for AmpliFinder tests."""

import os
import random
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from amplifinder.utils.file_utils import ensure_dir

# =============================================================================
# Assertion helpers
# =============================================================================


# =============================================================================
# Test configuration
# =============================================================================

FIXTURES_DIR = Path(__file__).parent / "fixtures"

# Root for real test data (MATLAB outputs, FASTQ, breseq, etc.)
TEST_DATA_ROOT = Path(
    os.environ.get("AMPLIFINDER_TEST_ROOT", "/zdata/user-data/rkishony/AmpliFinder_test")
)
MATLAB_OUTPUT = TEST_DATA_ROOT / "AmpliFinderWorkspace" / "output"

# Allow keeping outputs by setting AMPLIFINDER_KEEP_TEST_OUTPUT=true/1
KEEP_TEST_OUTPUT = True

# External data paths (from isolates.xlsx)
FASTQ_PATH = Path(
    "/zdata/user-data/idan/small_projects/breseq_to_amplification_from_dropbox/morbidostat_sra"
)
BRESEQ_PATH = Path(
    "/zdata/user-data/idan/small_projects/breseq_to_amplification_from_dropbox/Breseq1"
)


# =============================================================================
# Real-data integration fixtures (used across test suites)
# =============================================================================

# Persistent cache directory for test genomes (survives across test runs)
TEST_GENOMES_CACHE = Path(__file__).parent / "test_outputs" / "integration" / "genomesDB"


@pytest.fixture(scope="session")
def u00096_genome():
    """Load U00096 reference genome (session-scoped, cached across all tests).
    
    Downloads from NCBI only if not already cached. All integration tests
    share the same genome instance.
    """
    from amplifinder.steps import GetRefGenomeStep
    
    ref_path = ensure_dir(TEST_GENOMES_CACHE)
    
    step = GetRefGenomeStep(
        ref_name="U00096",
        ref_path=ref_path,
        ncbi=True,
    )
    return step.run()


@pytest.fixture
def isolate_srr25242877():
    """Test isolate SRR25242877 configuration."""
    return {
        "iso_name": "SRR25242877",
        "ref_name": "U00096",
        "fastq_path": FASTQ_PATH / "SRR25242877",
        "breseq_path": BRESEQ_PATH / "SRR25242877",
        "matlab_output": MATLAB_OUTPUT / "SRR25242877",
    }


@pytest.fixture
def isolate_srr25242906():
    """Test isolate SRR25242906 (ancestor)."""
    return {
        "iso_name": "SRR25242906",
        "ref_name": "U00096",
        "fastq_path": FASTQ_PATH / "SRR25242906",
        "breseq_path": BRESEQ_PATH / "SRR25242906",
        "matlab_output": MATLAB_OUTPUT / "SRR25242906",
    }


@pytest.fixture
def cleared_output_dir():
    """Clear entire output directory before test run.

    Returns the cleared output_dir Path.
    """
    base = Path(__file__).parent / "test_outputs" / "integration"
    output_dir = base / "output"

    ensure_dir(output_dir, cleanup=not KEEP_TEST_OUTPUT)

    return output_dir


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
            "topology": "circular",
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
    data_dir = ensure_dir(tmp_path / "test_data")

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
    out = ensure_dir(tmp_path / "output")
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


# =============================================================================
# Shared step fixtures
# =============================================================================

@pytest.fixture
def tiny_genome(tiny_ref_gbk, tiny_ref_fasta):
    """Create Genome object for tiny reference."""
    from amplifinder.data_types.genome import Genome
    return Genome(name="tiny", genbank_path=tiny_ref_gbk, fasta_path=tiny_ref_fasta)


@pytest.fixture
def locate_tns_step(tmp_output, tiny_genome):
    """Create TN location step (genbank-based)."""
    from amplifinder.steps import LocateTNsUsingGenbankStep
    return LocateTNsUsingGenbankStep(
        genome=tiny_genome,
        output_dir=tmp_output / "tn_loc" / tiny_genome.name,
    )


# =============================================================================
# Shared record fixtures
# =============================================================================

@pytest.fixture
def raw_tnjc2_record():
    """Base RawTnJc2 for step fixtures."""
    from amplifinder.data_types import RawTnJc2, Orientation

    return RawTnJc2(
        jc_num_S=1,
        jc_num_E=2,
        scaf="tiny",
        start=200,
        end=300,
        pos_tn_S=10,
        pos_tn_E=20,
        dir_tn_S=Orientation.FORWARD,
        dir_tn_E=Orientation.REVERSE,
        tn_ids=[1],
        tn_orientations=[Orientation.FORWARD],
        amplicon_length=100,
    )


@pytest.fixture
def covered_tnjc2_record(raw_tnjc2_record):
    """RawTnJc2 with coverage fields."""
    from amplifinder.data_types import CoveredTnJc2

    return CoveredTnJc2.from_other(
        raw_tnjc2_record,
        ref_name="tiny",
        iso_name="sample1",
        iso_scaf_avg=1.0,
        iso_amplicon_avg=2.0,
    )


@pytest.fixture
def classified_tnjc2_record(covered_tnjc2_record):
    """CoveredTnJc2 with structural classification."""
    from amplifinder.data_types import ClassifiedTnJc2, RawEvent

    return ClassifiedTnJc2.from_other(
        covered_tnjc2_record,
        raw_event=RawEvent.FLANKED,
        shared_tn_ids=[1],
        chosen_tn_id=1,
    )


@pytest.fixture
def filtered_tnjc2_record(classified_tnjc2_record):
    """FilteredTnJc2 with analysis directory."""
    from amplifinder.data_types import FilteredTnJc2

    return FilteredTnJc2.from_other(
        classified_tnjc2_record,
        analysis_dir="jc_200_300_001_L150",
    )


@pytest.fixture
def analyzed_tnjc2_record(filtered_tnjc2_record):
    """AnalyzedTnJc2 with junction coverage."""
    from amplifinder.data_types import AnalyzedTnJc2, RawEvent

    return AnalyzedTnJc2.from_other(
        filtered_tnjc2_record,
        jc_cov_left=[0] * 7,
        jc_cov_right=[0] * 7,
        jc_cov_spanning=[0] * 7,
        isolate_architecture=RawEvent.FLANKED,
        ancestor_architecture=None,
        event="flanked",
        event_modifiers=[],
    )


@pytest.fixture
def ref_tn_loc_record():
    """RefTnLoc base record."""
    from amplifinder.data_types import Orientation, RefTnLoc

    return RefTnLoc(
        tn_id=1,
        tn_name="IS1",
        tn_scaf="tiny",
        loc_left=200,
        loc_right=300,
        orientation=Orientation.FORWARD,
        join=False,
    )


def pytest_configure(config):
    """Register custom markers."""
    config.addinivalue_line(
        "markers", "integration: mark test as integration test (requires real data)"
    )
    config.addinivalue_line(
        "markers", "slow: mark test as slow (e.g., runs breseq)"
    )
