"""Pytest fixtures for AmpliFinder tests."""

import os
import random
import pandas as pd
import pytest

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from pathlib import Path

from amplifinder.data_types import RefTn, RecordTypedDf
from amplifinder.utils.file_utils import ensure_dir

# =============================================================================
# Test configuration
# =============================================================================

FIXTURES_DIR = Path(__file__).parent / "fixtures"

# Root for real test data (MATLAB outputs, FASTQ, breseq, etc.)
TEST_DATA_ROOT = Path(
    os.environ.get("AMPLIFINDER_TEST_ROOT", "/zdata/user-data/rkishony/AmpliFinder_test")
)
MATLAB_OUTPUT = TEST_DATA_ROOT / "AmpliFinderWorkspace" / "output"

# Output cleanup mode: "clear", "keep", or "keep_junctions"
OUTPUT_CLEANUP_MODE = "keep_junctions"

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


def _clear_except(directory: Path, exclude: list[str]):
    """Recursively delete all files/dirs except 'junctions' folders."""
    for item in directory.iterdir():
        if item.name in exclude:
            continue  # skip items in exclude
        if item.is_dir():
            _clear_except(item, exclude)
            # remove dir if empty after cleanup
            if not any(item.iterdir()):
                item.rmdir()
        else:
            item.unlink()


@pytest.fixture
def cleared_output_dir():
    """Clear output directory based on OUTPUT_CLEANUP_MODE.

    Modes:
        - "clear": Remove entire output directory
        - "keep": Keep all existing outputs
        - "keep_junctions": Clear all except junctions folder

    Returns the output_dir Path.
    """
    base = Path(__file__).parent / "test_outputs" / "integration"
    output_dir = base / "output"

    if OUTPUT_CLEANUP_MODE == "clear":
        ensure_dir(output_dir, cleanup=True)
    elif OUTPUT_CLEANUP_MODE == "keep":
        ensure_dir(output_dir, cleanup=False)
    elif OUTPUT_CLEANUP_MODE == "keep_junctions":
        if output_dir.exists():
            _clear_except(output_dir, exclude=['junctions'])
        else:
            ensure_dir(output_dir)
    else:
        raise ValueError(f"Invalid OUTPUT_CLEANUP_MODE: {OUTPUT_CLEANUP_MODE}")

    return output_dir


# =============================================================================
# Synthetic test data generators
# =============================================================================

def make_random_seq(length: int, seed: int = 42) -> Seq:
    """Generate reproducible random DNA sequence."""
    rng = random.Random(seed)
    return Seq("".join(rng.choices("ACGT", k=length)))


def make_tn_sequences() -> dict:
    """Define TN element sequences for testing.

    Uses real IS sequences from NCBI stored in fixtures directory.
    """
    is1_file = FIXTURES_DIR / "IS1.fasta"
    is5_file = FIXTURES_DIR / "IS5.fasta"

    is1_record = SeqIO.read(is1_file, "fasta")
    is5_record = SeqIO.read(is5_file, "fasta")

    return {
        "IS_test1": is1_record.seq,
        "IS_test2": is5_record.seq,
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
def locate_tns_step_factory(tmp_output, tiny_genome):
    """Factory for TN location steps (genbank-based)."""
    from amplifinder.steps import LocateTNsUsingGenbankStep

    def make():
        return LocateTNsUsingGenbankStep(
            genome=tiny_genome,
            output_dir=tmp_output / "tn_loc" / tiny_genome.name,
        )

    return make


@pytest.fixture
def locate_tns_step(locate_tns_step_factory):
    """Single TN location step instance (kept for compatibility)."""
    return locate_tns_step_factory()


# =============================================================================
# Shared record fixtures
# =============================================================================

@pytest.fixture
def raw_tnjc2_record(tiny_genome, ref_tn_record):
    """Base RawTnJc2 for step fixtures."""
    from amplifinder.data_types import (
        RawTnJc2,
        TnJunction,
        Orientation,
        OffsetRefTnSide,
        RefTnSide,
        Terminal,
        BaseEvent,
    )

    tn_jc_S = TnJunction(
        num=1,
        scaf1="tiny",
        pos1=10,
        dir1=Orientation.FORWARD,
        scaf2="tiny",
        pos2=200,
        dir2=Orientation.FORWARD,
        flanking1=50,
        flanking2=50,
        ref_tn_side=RefTnSide(ref_tn=ref_tn_record, side=Terminal.START),
        ref_tn_sides=[OffsetRefTnSide(ref_tn=ref_tn_record, side=Terminal.START, offset=0)],
        swapped=False,
    )
    tn_jc_E = TnJunction(
        num=2,
        scaf1="tiny",
        pos1=20,
        dir1=Orientation.REVERSE,
        scaf2="tiny",
        pos2=300,
        dir2=Orientation.REVERSE,
        flanking1=50,
        flanking2=50,
        ref_tn_side=RefTnSide(ref_tn=ref_tn_record, side=Terminal.END),
        ref_tn_sides=[OffsetRefTnSide(ref_tn=ref_tn_record, side=Terminal.END, offset=0)],
        swapped=False,
    )
    scaffold = tiny_genome.get_scaffold("tiny")
    return RawTnJc2(
        tnjc_left=tn_jc_S,
        tnjc_right=tn_jc_E,
        scaffold=scaffold,
        base_event=BaseEvent.REFERENCE_TN,
    )


@pytest.fixture
def covered_tnjc2_record(classified_tnjc2_record):
    """CoveredTnJc2 with coverage fields."""
    from amplifinder.data_types import CoveredTnJc2

    return CoveredTnJc2.from_other(
        classified_tnjc2_record,
        iso_scaf_avg=1.0,
        iso_amplicon_avg=2.0,
    )


@pytest.fixture
def classified_tnjc2_record(raw_tnjc2_record):
    """SingleLocusLinkedTnJc2 with structural classification."""
    from amplifinder.data_types import SingleLocusLinkedTnJc2

    return SingleLocusLinkedTnJc2.from_other(
        raw_tnjc2_record,
        single_locus_tnjc2_left_matchings=[],
        single_locus_tnjc2_right_matchings=[],
    )


@pytest.fixture
def covered_classified_tnjc2_record(classified_tnjc2_record):
    """CoveredTnJc2 with structural classification and coverage."""
    from amplifinder.data_types import CoveredTnJc2

    return CoveredTnJc2.from_other(
        classified_tnjc2_record,
        iso_scaf_avg=1.0,
        iso_amplicon_avg=2.0,
    )


@pytest.fixture
def filtered_tnjc2_record(covered_classified_tnjc2_record):
    """SynJctsTnJc2 with analysis directory."""
    from amplifinder.data_types import SynJctsTnJc2
    from amplifinder.data_types.rudimentary_junctions import RudimentaryJunctionValues
    from amplifinder.data_types.ref_tn import Terminal

    rudimentary = RudimentaryJunctionValues(
        amp_left_pos=200, amp_right_pos=300, amp_scaf="test_scaf",
        tn_start_pos=1, tn_end_pos=1000, tn_scaf="test_tn",
        tn_side_left_amp_side=Terminal.START,
        chr_left_pos_offset=0, chr_right_pos_offset=0,
        chr_left_pos_for_tn_offset=0, chr_right_pos_for_tn_offset=0,
        flank=150
    )

    return SynJctsTnJc2.from_other(
        covered_classified_tnjc2_record,
        rudimentary_junction_values=rudimentary,
    )


@pytest.fixture
def analyzed_tnjc2_record(filtered_tnjc2_record):
    """AnalyzedTnJc2 with junction coverage."""
    from amplifinder.data_types import AnalyzedTnJc2, JunctionType, JunctionReadCounts

    # Create jc_covs dict with all junction types having zero counts
    jc_covs = {jt: JunctionReadCounts() for jt in JunctionType}

    from amplifinder.data_types.rudimentary_junctions import RudimentaryJunctionValues
    from amplifinder.data_types.ref_tn import Terminal

    anc_rudimentary = RudimentaryJunctionValues(
        amp_left_pos=200, amp_right_pos=300, amp_scaf="test_scaf",
        tn_start_pos=1, tn_end_pos=1000, tn_scaf="test_tn",
        tn_side_left_amp_side=Terminal.START,
        chr_left_pos_offset=0, chr_right_pos_offset=0,
        chr_left_pos_for_tn_offset=0, chr_right_pos_for_tn_offset=0,
        flank=150
    )

    return AnalyzedTnJc2.from_other(
        filtered_tnjc2_record,
        anc_rudimentary_junction_values=anc_rudimentary,
        jc_covs=jc_covs,
        jc_covs_anc=None,
        jc_calls=None,
        jc_calls_anc=None,
    )


@pytest.fixture
def ref_tn_record():
    """RefTn base record."""
    from amplifinder.data_types import Orientation, RefTn

    return RefTn(
        tn_id=1,
        tn_name="IS1",
        scaf="tiny",
        is_circular=False,
        length=1000,
        start=200,
        end=300,
        orientation=Orientation.FORWARD,
        join=False,
    )


@pytest.fixture
def ref_tns_indexed(ref_tn_record):
    """RefTn RecordTypedDf with tn_id as index (for O(1) lookup in synthetic junctions)."""

    ref_tn_df = pd.DataFrame([ref_tn_record.model_dump()])
    ref_tn_df = ref_tn_df.set_index('tn_id', drop=False)  # Keep tn_id column
    return RecordTypedDf(ref_tn_df, RefTn)


def pytest_configure(config):
    """Register custom markers."""
    config.addinivalue_line(
        "markers", "integration: mark test as integration test (requires real data)"
    )
    config.addinivalue_line(
        "markers", "slow: mark test as slow (e.g., runs breseq)"
    )


@pytest.fixture(autouse=True)
def disable_rich_console_for_tests(monkeypatch):
    """Disable rich console colors for testing so caplog/capsys work properly."""
    from amplifinder.logger import logger
    monkeypatch.setattr(logger, 'use_colors', False)
    monkeypatch.setattr(logger, 'console', None)
