"""Integration tests comparing Python pipeline to MATLAB output."""

import pandas as pd
import pytest

from amplifinder.steps import (
    GetReferenceStep,
    LocateTNsUsingGenbankStep,
    LocateTNsUsingISfinderStep,
    CreateReferenceTnJunctionsStep,
    CreateRefTnEndSeqsStep,
    CreateTNJCStep,
    CreateTNJC2Step,
)
from amplifinder.data_types import RecordTypedDF, Junction

pytestmark = pytest.mark.integration


# =============================================================================
# Test: Reference loading
# =============================================================================

class TestGetReference:
    """Test reference genome loading for U00096."""

    def test_load_u00096_from_ncbi(self, tmp_path, isolate_srr25242877):
        """Load E. coli U00096 reference from NCBI."""
        ref_path = tmp_path / "genomesDB"
        ref_path.mkdir()

        step = GetReferenceStep(
            ref_name="U00096",
            ref_path=ref_path,
            ncbi=True,
        )
        genome = step.run()

        assert genome.name == "U00096"
        assert genome.genbank_path.exists()
        assert genome.fasta_path.exists()
        assert genome.length > 4_000_000  # E. coli ~4.6M bp


# =============================================================================
# Test: TN location
# =============================================================================

class TestLocateTNs:
    """Test TN element location using GenBank annotations."""

    @pytest.fixture
    def u00096_genome(self, tmp_path):
        """Load U00096 reference (cached)."""
        ref_path = tmp_path / "genomesDB"
        ref_path.mkdir()

        step = GetReferenceStep(
            ref_name="U00096",
            ref_path=ref_path,
            ncbi=True,
        )
        return step.run()

    def test_locate_tns_genbank(self, tmp_path, u00096_genome, isolate_srr25242877):
        """Locate TN elements from GenBank annotations."""
        step = LocateTNsUsingGenbankStep(
            genbank_path=u00096_genome.genbank_path,
            ref_name=u00096_genome.name,
            ref_path=tmp_path,
        )
        tn_loc = step.run()

        # MATLAB found IS elements in U00096
        assert tn_loc is not None
        assert len(tn_loc) > 0

        # Check expected columns
        assert "TN_Name" in tn_loc.df.columns
        assert "LocLeft" in tn_loc.df.columns
        assert "LocRight" in tn_loc.df.columns

    def test_locate_tns_isfinder(self, tmp_path, u00096_genome):
        """Locate TN elements using ISfinder database."""
        from amplifinder.data import get_builtin_isfinder_db_path

        step = LocateTNsUsingISfinderStep(
            ref_fasta=u00096_genome.fasta_path,
            ref_name=u00096_genome.name,
            ref_path=tmp_path,
            isdb_path=get_builtin_isfinder_db_path(),
            evalue=1e-4,
            critical_coverage=0.9,
        )
        tn_loc = step.run()

        assert len(tn_loc) > 0


# =============================================================================
# Test: Breseq parsing (uses pre-computed breseq output)
# =============================================================================

class TestBreseq:
    """Test breseq output parsing using existing breseq results."""

    def test_parse_existing_breseq_output(self, tmp_path, isolate_srr25242877):
        """Parse pre-computed breseq output for SRR25242877."""
        breseq_path = isolate_srr25242877["breseq_path"]

        if not breseq_path.exists():
            pytest.skip(f"Breseq output not found: {breseq_path}")

        # Parse the output.gd file directly
        from amplifinder.tools.breseq import parse_breseq_output

        output_gd = breseq_path / "output" / "output.gd"
        if not output_gd.exists():
            pytest.skip(f"output.gd not found: {output_gd}")

        # parse_breseq_output expects the breseq output directory, not the .gd file
        result = parse_breseq_output(breseq_path)

        # Check JC (junction) output
        assert "JC" in result
        jc_df = result["JC"]
        assert len(jc_df) > 0


# =============================================================================
# Test: Compare with MATLAB output
# =============================================================================

class TestCompareWithMATLAB:
    """Compare Python outputs with MATLAB reference outputs."""

    def test_isjc2_structure(self, isolate_srr25242877):
        """Verify we can load and parse MATLAB ISJC2.xlsx output."""
        matlab_path = isolate_srr25242877["matlab_output"]
        isjc2_xlsx = matlab_path / "ISJC2.xlsx"

        if not isjc2_xlsx.exists():
            pytest.skip(f"MATLAB output not found: {isjc2_xlsx}")

        df = pd.read_excel(isjc2_xlsx)

        # Expected columns from MATLAB output
        expected_cols = [
            "iso", "Reference", "IS_element", "IS_direction",
            "amplicon_length", "median_copy_number",
        ]
        for col in expected_cols:
            assert col in df.columns, f"Missing column: {col}"

    def test_classified_amplifications(self, isolate_srr25242877):
        """Verify MATLAB classified_amplifications output."""
        matlab_path = isolate_srr25242877["matlab_output"]
        xlsx = matlab_path / "classified_amplifications.xlsx"

        if not xlsx.exists():
            pytest.skip(f"MATLAB output not found: {xlsx}")

        df = pd.read_excel(xlsx)
        assert len(df) >= 1  # At least one amplification found


# =============================================================================
# Full pipeline comparison (slow)
# =============================================================================

@pytest.mark.slow
class TestFullPipeline:
    """Full pipeline integration tests."""

    def test_pipeline_matches_matlab_tn_count(self, tmp_path, isolate_srr25242877):
        """Run pipeline and compare TN junction count with MATLAB."""
        matlab_path = isolate_srr25242877["matlab_output"]
        isjc2_xlsx = matlab_path / "ISJC2.xlsx"

        if not isjc2_xlsx.exists():
            pytest.skip("MATLAB reference output not available")

        matlab_df = pd.read_excel(isjc2_xlsx)
        matlab_jc_count = len(matlab_df)

        # TODO: Run full pipeline and compare
        # For now, just verify MATLAB output is accessible
        assert matlab_jc_count > 0
        print(f"MATLAB found {matlab_jc_count} TN junctions")

    def test_tnjc2_step_produces_output(self, tmp_path, isolate_srr25242877):
        """Test TNJC2 step produces valid output structure."""
        breseq_path = isolate_srr25242877["breseq_path"]

        if not breseq_path.exists():
            pytest.skip(f"Breseq output not found: {breseq_path}")

        ref_path = tmp_path / "genomesDB"
        ref_path.mkdir()
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        # Step 1: Get reference
        genome = GetReferenceStep(
            ref_name="U00096",
            ref_path=ref_path,
            ncbi=True,
        ).run()

        # Step 2: Locate TNs
        tn_loc = LocateTNsUsingGenbankStep(
            genbank_path=genome.genbank_path,
            ref_name=genome.name,
            ref_path=ref_path,
        ).run()

        assert tn_loc is not None and len(tn_loc) > 0

        # Step 3: Create reference TN junctions
        ref_tn_jc = CreateReferenceTnJunctionsStep(
            tn_loc=tn_loc,
            ref_name=genome.name,
            output_dir=output_dir,
            reference_tn_out_span=50,
        ).run()

        ref_tn_end_seqs = CreateRefTnEndSeqsStep(
            ref_tn_jc=ref_tn_jc,
            genome=genome,
            ref_path=ref_path,
            max_dist_to_tn=200,
        ).run()

        # Step 4: Parse breseq output
        from amplifinder.tools.breseq import parse_breseq_output
        breseq_output = parse_breseq_output(breseq_path)
        breseq_jc = breseq_output["JC"]

        # Step 5: Create TNJC
        all_jc_df = pd.concat([breseq_jc, ref_tn_jc.df], ignore_index=True)
        all_jc = RecordTypedDF(all_jc_df, Junction)

        tnjc = CreateTNJCStep(
            jc_df=all_jc,
            ref_tn_end_seqs=ref_tn_end_seqs,
            genome=genome,
            output_dir=output_dir,
            max_dist_to_tn=200,
            trim_jc_flanking=5,
        ).run()

        # Step 6: Create TNJC2
        tnjc2 = CreateTNJC2Step(
            tnjc=tnjc,
            genome=genome,
            output_dir=output_dir,
        ).run()

        # Verify output structure
        assert isinstance(tnjc2, RecordTypedDF)
        assert (output_dir / "TNJC2.csv").exists()

        # Check expected columns
        expected_cols = {
            "jc_num_L", "jc_num_R", "scaf_chr",
            "amplicon_length", "complementary_length",
            "tn_orientation", "span_origin",
        }
        assert expected_cols.issubset(set(tnjc2.df.columns))

        print(f"TNJC2: found {len(tnjc2)} junction pairs")
