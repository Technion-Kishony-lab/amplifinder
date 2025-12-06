"""Tests for BLAST tools."""

import pytest

from amplifinder.tools.blast import run_blastn, parse_blast_csv, make_blast_db
from amplifinder.data.schemas import BLAST_SCHEMA
from tests.env import RUN_BLAST_TESTS


skip_no_blast = pytest.mark.skipif(not RUN_BLAST_TESTS, reason="BLAST tests disabled")


@pytest.fixture
def blast_db(tiny_tn_db, tmp_path):
    """Create BLAST database and return path."""
    db_path = tmp_path / "test_db"
    make_blast_db(fasta=tiny_tn_db, db_path=db_path)
    return db_path


@pytest.fixture
def blast_output(tiny_ref_fasta, blast_db, tmp_path):
    """Run BLAST and return output path."""
    output = tmp_path / "blast_out.txt"
    run_blastn(query=tiny_ref_fasta, db=blast_db, output=output)
    return output


@skip_no_blast
def test_make_blast_db(blast_db):
    """Should create BLAST database files."""
    db_dir = blast_db.parent
    assert (db_dir / "test_db.nhr").exists() or (db_dir / "test_db.nsq").exists()


@skip_no_blast
def test_run_blastn(blast_output):
    """Should run BLAST and produce output with hits."""
    assert blast_output.exists()
    assert len(blast_output.read_text()) > 0


@skip_no_blast
def test_parse_blast_csv(blast_output):
    """Should parse BLAST CSV output correctly."""
    df = parse_blast_csv(blast_output)

    assert len(df) > 0
    BLAST_SCHEMA.assert_matches(df)


@skip_no_blast
def test_blast_no_hits(blast_db, tmp_path):
    """Should return empty DataFrame with correct schema when no hits."""
    # Create query with random sequence that won't match TN database
    from tests.conftest import make_random_seq, write_fasta
    query = tmp_path / "no_match.fasta"
    write_fasta(query, "no_match", make_random_seq(200, seed=999))

    output = tmp_path / "blast_out.txt"
    run_blastn(query=query, db=blast_db, output=output)
    df = parse_blast_csv(output)

    assert len(df) == 0
    BLAST_SCHEMA.assert_matches(df)
