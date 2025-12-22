"""Integration tests comparing Python pipeline to MATLAB output."""

import pandas as pd
import pytest

from amplifinder.steps import (
    GetRefGenomeStep,
    LocateTNsUsingGenbankStep,
    LocateTNsUsingISfinderStep,
    CreateRefTnJcStep,
    CreateRefTnEndSeqsStep,
    CreateTnJcStep,
    CreateTnJc2Step,
)
from amplifinder.data_types import RecordTypedDF, Junction
from amplifinder.pipeline import Pipeline

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

        step = GetRefGenomeStep(
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

        step = GetRefGenomeStep(
            ref_name="U00096",
            ref_path=ref_path,
            ncbi=True,
        )
        return step.run()

    def test_locate_tns_genbank(self, tmp_path, u00096_genome, isolate_srr25242877):
        """Locate TN elements from GenBank annotations."""
        output_dir = tmp_path / "tn_loc" / u00096_genome.name
        step = LocateTNsUsingGenbankStep(
            genome=u00096_genome,
            output_dir=output_dir,
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

        output_dir = tmp_path / "tn_loc" / u00096_genome.name
        step = LocateTNsUsingISfinderStep(
            genome=u00096_genome,
            output_dir=output_dir,
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
        # Use tmp_path for CSV output to avoid permission errors
        result = parse_breseq_output(breseq_path, csv_output_dir=tmp_path)

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

    def test_full_pipeline_matches_matlab(self, isolate_srr25242877):
        """Run full pipeline and compare with MATLAB outputs (1-to-1 matching)."""
        from amplifinder.config import Config, get_run_dir
        from amplifinder.pipeline import Pipeline
        from tests.test_integration.matlab_compare import compare_isjc2_outputs
        
        # Setup output directory next to MATLAB outputs
        matlab_output_dir = isolate_srr25242877["matlab_output"]
        # matlab_output_dir = /zdata/user-data/rkishony/AmpliFinder_test/AmpliFinderWorkspace/output/SRR25242877
        # Go up to AmpliFinder_test, then create python_outputs directory
        test_output_root = matlab_output_dir.parent.parent.parent / "python_outputs"
        # test_output_root = /zdata/user-data/rkishony/AmpliFinder_test/python_outputs
        
        # Check MATLAB output exists
        isjc2_xlsx = matlab_output_dir / "ISJC2.xlsx"
        if not isjc2_xlsx.exists():
            pytest.skip("MATLAB reference output not available")
        
        # Setup config using pre-computed breseq output
        config = Config(
            iso_path=isolate_srr25242877["fastq_path"],
            ref_name="U00096",
            anc_path=None,  # Can add ancestor later
            iso_name="SRR25242877",
            output_dir=test_output_root / "output",
            ref_path=test_output_root / "genomesDB",
            iso_breseq_path=isolate_srr25242877["breseq_path"],
            ncbi=True,
            use_isfinder=False,  # Use GenBank for consistency
        )
        
        # Run full pipeline
        pipeline = Pipeline(config)
        result = pipeline.run()
        
        # Verify outputs exist
        run_dir = get_run_dir(config)
        # run_dir = test_output_root / "output" / "U00096" / "SRR25242877" / "SRR25242877"
        
        # Check if pipeline produced any candidates
        if len(result) == 0:
            pytest.fail(
                f"Pipeline found 0 candidates. Expected to match MATLAB output with 155 junctions. "
                f"Check pipeline parameters and matching logic. Run dir: {run_dir}"
            )
        
        assert (run_dir / "ISJC2.csv").exists(), f"ISJC2.csv not found in {run_dir}"
        assert (run_dir / "candidate_amplifications.csv").exists(), f"candidate_amplifications.csv not found in {run_dir}"
        
        # Load Python outputs
        python_isjc2 = pd.read_csv(run_dir / "ISJC2.csv")
        print(f"Python found {len(python_isjc2)} candidates")
        
        # Load MATLAB outputs
        matlab_isjc2 = pd.read_excel(isjc2_xlsx)
        print(f"MATLAB found {len(matlab_isjc2)} candidates")
        
        # Compare outputs (1-to-1 matching required)
        compare_isjc2_outputs(python_isjc2, matlab_isjc2)

    def test_tnjc2_step_produces_output(self, tmp_path, isolate_srr25242877):
        """Test TnJc2 step produces valid output structure."""
        breseq_path = isolate_srr25242877["breseq_path"]

        if not breseq_path.exists():
            pytest.skip(f"Breseq output not found: {breseq_path}")

        ref_path = tmp_path / "genomesDB"
        ref_path.mkdir()
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        # Step 1: Get reference
        genome = GetRefGenomeStep(
            ref_name="U00096",
            ref_path=ref_path,
            ncbi=True,
        ).run()

        # Step 2: Locate TNs
        tn_loc = LocateTNsUsingGenbankStep(
            genome=genome,
            output_dir=ref_path / "tn_loc" / genome.name,
        ).run()

        assert tn_loc is not None and len(tn_loc) > 0

        # Step 3: Create reference TN junctions
        ref_tn_jc = CreateRefTnJcStep(
            tn_loc=tn_loc,
            genome=genome,
            output_dir=output_dir,
            source="isfinder",
            reference_tn_out_span=50,
        ).run()

        ref_tn_end_seqs = CreateRefTnEndSeqsStep(
            ref_tn_jc=ref_tn_jc,
            genome=genome,
            output_dir=ref_path / "tn_loc" / genome.name,
            source="isfinder",
            max_dist_to_tn=200,
        ).run()

        # Step 4: Parse breseq output
        from amplifinder.tools.breseq import parse_breseq_output
        # Use tmp_path for CSV output to avoid permission errors
        breseq_output = parse_breseq_output(breseq_path, csv_output_dir=tmp_path)
        breseq_jc = breseq_output["JC"]

        # Step 5: Create TnJc
        all_jc_df = pd.concat([breseq_jc, ref_tn_jc.df], ignore_index=True)
        all_jc = RecordTypedDF(all_jc_df, Junction)

        tnjc = CreateTnJcStep(
            jc_df=all_jc,
            ref_tn_end_seqs=ref_tn_end_seqs,
            genome=genome,
            output_dir=output_dir,
            max_dist_to_tn=200,
            trim_jc_flanking=5,
        ).run()

        # Step 6: Create TnJc2
        tnjc2 = CreateTnJc2Step(
            tnjc=tnjc,
            genome=genome,
            output_dir=output_dir,
        ).run()

        # Verify output structure
        assert isinstance(tnjc2, RecordTypedDF)
        assert (output_dir / "tn_jc2.csv").exists()

        # Check expected columns
        expected_cols = {
            "jc_num_L", "jc_num_R", "scaf_chr",
            "amplicon_length", "complementary_length",
            "tn_orientations", "span_origin",
        }
        assert expected_cols.issubset(set(tnjc2.df.columns))

        print(f"TnJc2: found {len(tnjc2)} junction pairs")


# =============================================================================
# Step-by-step debugging tests (slow) - DIAGNOSTIC TOOL
# =============================================================================

class TestPipeline(Pipeline):
    """Pipeline subclass for testing with step-by-step MATLAB comparisons."""
    
    def __init__(self, config, matlab_output_dir):
        super().__init__(config)
        self.matlab_output_dir = matlab_output_dir
        self.step_results = {}
        self.is_ancestor_run = False
    
    def _locate_tns_in_reference(self, genome):
        """Step 2: Locate TN elements - compare with MATLAB."""
        result = super()._locate_tns_in_reference(genome)
        self.step_results["tn_loc"] = result
        
        # Compare with MATLAB
        import scipy.io as sio
        matlab_file = self.matlab_output_dir / "genbank_IS_loc.mat"
        if matlab_file.exists():
            matlab_data = sio.loadmat(str(matlab_file))
            matlab_is_loc = matlab_data.get("IS_loc", None)
            matlab_count = len(matlab_is_loc) if matlab_is_loc is not None else 0
            print(f"Step 2: MATLAB={matlab_count} TNs, Python={len(result)} TNs")
        else:
            print(f"Step 2: Python={len(result)} TNs (MATLAB file not found: {matlab_file})")
        
        return result
    
    def _create_reference_tn_junctions(self, tn_loc, genome, iso_output):
        """Step 3: Create reference TN junctions - compare with MATLAB."""
        result = super()._create_reference_tn_junctions(tn_loc, genome, iso_output)
        ref_tn_jc, ref_tn_end_seqs = result
        self.step_results["ref_tn_jc"] = ref_tn_jc
        self.step_results["ref_tn_end_seqs"] = ref_tn_end_seqs
        
        # Compare with MATLAB
        import scipy.io as sio
        matlab_file = self.matlab_output_dir / "JC_refIS.mat"
        if matlab_file.exists():
            matlab_data = sio.loadmat(str(matlab_file))
            matlab_jc_refis = matlab_data.get("JC_refIS", None)
            matlab_count = len(matlab_jc_refis) if matlab_jc_refis is not None else 0
            print(f"Step 3: MATLAB={matlab_count} ref junctions, Python={len(ref_tn_jc)} ref junctions")
        else:
            print(f"Step 3: Python={len(ref_tn_jc)} ref junctions (MATLAB file not found: {matlab_file})")
        
        return result
    
    def _run_breseq(self, genome, iso_output):
        """Step 4: Parse breseq - compare with MATLAB."""
        result = super()._run_breseq(genome, iso_output)
        self.step_results["breseq_jc"] = result
        
        # Compare with MATLAB
        import scipy.io as sio
        matlab_file = self.matlab_output_dir / "MUT.mat"
        if matlab_file.exists():
            matlab_data = sio.loadmat(str(matlab_file))
            matlab_jc = matlab_data.get("JC", None)
            matlab_count = len(matlab_jc) if matlab_jc is not None else 0
            print(f"Step 4: MATLAB={matlab_count} junctions, Python={len(result)} junctions")
        else:
            print(f"Step 4: Python={len(result)} junctions (MATLAB file not found: {matlab_file})")
        
        return result
    
    def _create_tnjc(self, breseq_jc, ref_tn_jc, ref_tn_end_seqs, genome, iso_output):
        """Step 5: Match junctions to TN elements - compare with MATLAB (WHERE IT FAILS)."""
        result = super()._create_tnjc(breseq_jc, ref_tn_jc, ref_tn_end_seqs, genome, iso_output)
        self.step_results["tnjc"] = result
        
        # Compare with MATLAB
        import scipy.io as sio
        matlab_file = self.matlab_output_dir / "ISJC.mat"
        if matlab_file.exists():
            matlab_data = sio.loadmat(str(matlab_file))
            matlab_isjc = matlab_data.get("ISJC", None)
            matlab_count = len(matlab_isjc) if matlab_isjc is not None else 0
            print(f"Step 5: MATLAB={matlab_count} IS-associated junctions, Python={len(result)} TN-associated junctions")
            
            if len(result) == 0 and matlab_count > 0:
                print("\n*** DEBUGGING NEEDED ***")
                print("MATLAB found IS-associated junctions but Python found 0")
                print("Check:")
                print("  1. Junction sequences extraction")
                print("  2. TN end sequences")
                print("  3. Matching parameters (max_dist_to_IS)")
                print("  4. Scaffold name matching")
        else:
            print(f"Step 5: Python={len(result)} TN-associated junctions (MATLAB file not found: {matlab_file})")
        
        return result
    
    def _create_tnjc2(self, tnjc, genome, iso_output):
        """Step 6: Combine junction pairs - compare with MATLAB."""
        result = super()._create_tnjc2(tnjc, genome, iso_output)
        self.step_results["tnjc2"] = result
        
        # Compare with MATLAB
        import scipy.io as sio
        matlab_file = self.matlab_output_dir / "ISJC2.mat"
        if matlab_file.exists():
            matlab_data = sio.loadmat(str(matlab_file))
            matlab_isjc2 = matlab_data.get("ISJC2", None)
            matlab_count = len(matlab_isjc2) if matlab_isjc2 is not None else 0
            print(f"Step 6: MATLAB={matlab_count} junction pairs, Python={len(result)} junction pairs")
        else:
            print(f"Step 6: Python={len(result)} junction pairs (MATLAB file not found: {matlab_file})")
        
        return result


@pytest.mark.slow
class TestPipelineStepByStep:
    """Step-by-step pipeline tests for debugging - compares each step with MATLAB."""
    
    def test_isolate_pipeline_steps_with_matlab(self, isolate_srr25242877, tmp_path):
        """Test isolate pipeline step-by-step, comparing each step with MATLAB outputs."""
        from amplifinder.config import Config
        
        matlab_output_dir = isolate_srr25242877["matlab_output"]
        test_output_root = matlab_output_dir.parent.parent.parent / "python_outputs"
        
        config = Config(
            iso_path=isolate_srr25242877["fastq_path"],
            ref_name="U00096",
            anc_path=None,  # No ancestor for this test
            iso_name="SRR25242877",
            output_dir=test_output_root / "output",
            ref_path=test_output_root / "genomesDB",
            iso_breseq_path=isolate_srr25242877["breseq_path"],
            ncbi=True,
            use_isfinder=False,  # Use GenBank for consistency
        )
        
        print("\n=== Testing Isolate Pipeline (no ancestor) ===")
        pipeline = TestPipeline(config, matlab_output_dir)
        result = pipeline.run()
        
        print(f"\nFinal result: {len(result)} candidates")
        return result
    
    def test_ancestor_pipeline_steps_with_matlab(self, isolate_srr25242906, tmp_path):
        """Test ancestor pipeline step-by-step, comparing each step with MATLAB outputs."""
        from amplifinder.config import Config
        
        matlab_output_dir = isolate_srr25242906["matlab_output"]
        test_output_root = matlab_output_dir.parent.parent.parent / "python_outputs"
        
        config = Config(
            iso_path=isolate_srr25242906["fastq_path"],
            ref_name="U00096",
            anc_path=None,  # Ancestor has no ancestor
            iso_name="SRR25242906",
            anc_name="SRR25242906",
            output_dir=test_output_root / "output",
            ref_path=test_output_root / "genomesDB",
            iso_breseq_path=isolate_srr25242906["breseq_path"],
            ncbi=True,
            use_isfinder=False,
        )
        
        print("\n=== Testing Ancestor Pipeline ===")
        pipeline = TestPipeline(config, matlab_output_dir)
        pipeline.is_ancestor_run = True
        result = pipeline.run()
        
        print(f"\nFinal result: {len(result)} candidates")
        return result
    
    def test_isolate_with_ancestor_pipeline_steps(self, isolate_srr25242877, isolate_srr25242906, tmp_path):
        """Test isolate pipeline with ancestor - step-by-step comparison."""
        from amplifinder.config import Config
        
        matlab_output_dir = isolate_srr25242877["matlab_output"]
        test_output_root = matlab_output_dir.parent.parent.parent / "python_outputs"
        
        config = Config(
            iso_path=isolate_srr25242877["fastq_path"],
            ref_name="U00096",
            anc_path=isolate_srr25242906["fastq_path"],  # Has ancestor
            iso_name="SRR25242877",
            anc_name="SRR25242906",
            anc_breseq_path=isolate_srr25242906["breseq_path"],
            output_dir=test_output_root / "output",
            ref_path=test_output_root / "genomesDB",
            iso_breseq_path=isolate_srr25242877["breseq_path"],
            ncbi=True,
            use_isfinder=False,
        )
        
        print("\n=== Testing Isolate Pipeline (with ancestor) ===")
        print("Note: Ancestor pipeline will run automatically if needed")
        pipeline = TestPipeline(config, matlab_output_dir)
        result = pipeline.run()  # Will run ancestor automatically if needed
        
        print(f"\nFinal result: {len(result)} candidates")
        return result
