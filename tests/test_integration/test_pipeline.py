"""Integration tests for full pipeline execution."""

import os
from pathlib import Path

import pandas as pd
import pytest

from amplifinder.steps import (
    GetRefGenomeStep,
    LocateTNsUsingGenbankStep,
    CreateRefTnJcStep,
    CreateRefTnEndSeqsStep,
    CreateTnJcStep,
    CreateTnJc2Step,
)
from amplifinder.steps.base import Step
from amplifinder.config import Config
from amplifinder.data_types import RecordTypedDf, Junction
from amplifinder.pipeline import Pipeline

pytestmark = pytest.mark.integration

# Global flag: if True, skip tests when MATLAB files are missing; if False, fail
REQUIRE_MATLAB_FILES = True


# =============================================================================
# Full pipeline comparison (slow)
# =============================================================================

@pytest.mark.slow
class TestPipelineStepByStep:
    """Pipeline integration tests with step-by-step MATLAB comparison."""


# =============================================================================
# Step-by-step debugging tests (slow) - DIAGNOSTIC TOOL
# =============================================================================

class TestPipeline(Pipeline):
    """Pipeline subclass for testing with step-by-step MATLAB comparisons."""
    
    def __init__(self, config, matlab_output_dir):
        super().__init__(config)
        self.matlab_output_dir = matlab_output_dir
    
    def _load_matlab_isjc2(self):
        """Load MATLAB ISJC2.xlsx (final - junction pairs) if available.
        
        If REQUIRE_MATLAB_FILES is True, fails if file doesn't exist.
        Otherwise returns None if file doesn't exist.
        """
        import pandas as pd
        isjc2_xlsx = self.matlab_output_dir / "ISJC2.xlsx"
        if isjc2_xlsx.exists():
            return pd.read_excel(isjc2_xlsx)
        
        if REQUIRE_MATLAB_FILES:
            pytest.fail(f"MATLAB ISJC2.xlsx not found: {isjc2_xlsx}")
        return None
    
    def _locate_tns_in_reference(self, genome):
        """Step 2: Locate TN elements - compare with MATLAB."""
        result = super()._locate_tns_in_reference(genome)
        print(f"Step 2: Python={len(result)} TNs")
        return result
    
    def _create_reference_tn_junctions(self, tn_loc, genome, iso_output):
        """Step 3: Create reference TN junctions - compare with MATLAB."""
        result = super()._create_reference_tn_junctions(tn_loc, genome, iso_output)
        ref_tn_jc, ref_tn_end_seqs = result
        print(f"Step 3: Python={len(ref_tn_jc)} ref junctions")
        return result
    
    def _run_breseq(self, genome, iso_output):
        """Step 4: Parse breseq - compare with MATLAB."""
        result = super()._run_breseq(genome, iso_output)
        print(f"Step 4: Python={len(result)} junctions")
        return result
    
    def _create_tnjc(self, breseq_jc, ref_tn_jc, ref_tn_end_seqs, genome, iso_output):
        """Step 5: Match junctions to TN elements."""
        result = super()._create_tnjc(breseq_jc, ref_tn_jc, ref_tn_end_seqs, genome, iso_output)
        print(f"Step 5: Python={len(result)} TN-associated junctions")
        return result
    
    def _create_tnjc2(self, tnjc, genome, iso_output):
        """Step 6: Combine junction pairs."""
        result = super()._create_tnjc2(tnjc, genome, iso_output)
        print(f"Step 6: Python={len(result)} junction pairs")
        return result
    
    def _calc_coverage(self, tnjc2, genome, iso_output):
        """Step 7: Calculate amplicon coverage - compare CoveredTnJc2 with MATLAB ISJC2."""
        result = super()._calc_amplicon_coverage(tnjc2, genome, iso_output)
        
        # Compare CoveredTnJc2 with MATLAB ISJC2.xlsx
        matlab_df = self._load_matlab_isjc2()
        if matlab_df is not None:
            print(f"Step 7: MATLAB ISJC2.xlsx has {len(matlab_df)} candidates, Python={len(result)} candidates")
            # Convert CoveredTnJc2 to comparable format
            python_df = result.df.copy()
            python_df['Positions_in_chromosome'] = python_df.apply(
                lambda row: f"{row['pos_chr_L']}-{row['pos_chr_R']}", axis=1
            )
            python_df['amplicon_length'] = python_df['amplicon_length']
            python_df['IS_element'] = python_df['tn_ids'].apply(
                lambda x: ','.join(map(str, x)) if isinstance(x, list) else str(x)
            )
            python_df['amplicon_coverage'] = python_df['amplicon_coverage']
            
            # Compare with MATLAB ISJC2
            from tests.test_integration.matlab_compare import compare_isjc2_outputs
            compare_isjc2_outputs(python_df, matlab_df)
            print(f"Step 7: ✓ Comparison passed: CoveredTnJc2 matches MATLAB ISJC2")
        else:
            assert not REQUIRE_MATLAB_FILES, "MATLAB file missing but REQUIRE_MATLAB_FILES=True"
            print(f"Step 7: Python={len(result)} candidates (MATLAB file not found)")
        
        return result
    
    def _export(self, analyzed, iso_output):
        """Step 14: Export results - compare final export tables with MATLAB."""
        super()._export(analyzed, iso_output)
        
        # Compare final export tables with MATLAB
        from tests.test_integration.matlab_compare import (
            load_matlab_candidate, load_matlab_classified, compare_isjc2_outputs
        )
        
        # Compare candidate_amplifications.csv with MATLAB candidate_amplifications.xlsx
        python_candidate_file = iso_output / "candidate_amplifications.csv"
        if python_candidate_file.exists():
            matlab_candidate = load_matlab_candidate(self.matlab_output_dir)
            if matlab_candidate is not None:
                python_candidate = pd.read_csv(python_candidate_file)
                print(f"\n=== Comparing candidate_amplifications ===")
                print(f"Python: {len(python_candidate)} candidates")
                print(f"MATLAB: {len(matlab_candidate)} candidates")
                compare_isjc2_outputs(python_candidate, matlab_candidate)
                print("✓ candidate_amplifications comparison passed")
            else:
                print(f"MATLAB candidate_amplifications.xlsx not found (skipping comparison)")
        
        # Compare classified_amplifications (if exists)
        # Note: Python creates ISJC2.csv which contains all analyzed candidates
        # MATLAB creates classified_amplifications.xlsx
        python_isjc2_file = iso_output / "ISJC2.csv"
        if python_isjc2_file.exists():
            matlab_classified = load_matlab_classified(self.matlab_output_dir)
            if matlab_classified is not None:
                python_isjc2 = pd.read_csv(python_isjc2_file)
                print(f"\n=== Comparing classified_amplifications ===")
                print(f"Python ISJC2.csv: {len(python_isjc2)} candidates")
                print(f"MATLAB classified_amplifications.xlsx: {len(matlab_classified)} candidates")
                compare_isjc2_outputs(python_isjc2, matlab_classified)
                print("✓ classified_amplifications comparison passed")
            else:
                print(f"MATLAB classified_amplifications.xlsx not found (skipping comparison)")


@pytest.mark.slow
class TestPipelineStepByStep:
    """Pipeline integration tests with step-by-step MATLAB comparison."""
    
    @staticmethod
    def _setup_pipeline(config, return_run_dir=False):
        """Common setup: enable verbose reporting."""
        Step.set_global_verbose(True)
        if return_run_dir:
            from amplifinder.config import get_iso_run_dir
            return get_iso_run_dir(config)
        return None
    
    @staticmethod
    def _get_test_output_root(matlab_output_dir):
        """Get test output root directory next to MATLAB outputs.

        By default this is:
            {AMPLIFINDER_TEST_ROOT or default}/python_outputs

        Override base directory with AMPLIFINDER_OUTPUT_ROOT.
        """
        default_base = matlab_output_dir.parent.parent.parent
        base = Path(os.environ.get("AMPLIFINDER_OUTPUT_ROOT", default_base))
        return base / "python_outputs"
    
    @staticmethod
    def _create_config(isolate, output_dir, anc_isolate=None, anc_name=None):
        """Create Config with common defaults."""
        test_output_root = output_dir.parent
        
        config_kwargs = {
            "iso_path": isolate["fastq_path"],
            "ref_name": "U00096",
            "iso_name": isolate["iso_name"],
            "output_dir": output_dir,
            "ref_path": test_output_root / "genomesDB",
            "iso_breseq_path": isolate["breseq_path"],
            "ncbi": True,
            "use_isfinder": False,
        }
        
        if anc_isolate is not None:
            config_kwargs["anc_path"] = anc_isolate["fastq_path"]
            config_kwargs["anc_name"] = anc_name or anc_isolate["iso_name"]
            config_kwargs["anc_breseq_path"] = anc_isolate["breseq_path"]
        else:
            config_kwargs["anc_path"] = None
        
        return Config(**config_kwargs)
    
    def test_full_pipeline_matches_matlab(self, isolate_srr25242877, cleared_output_dir):
        """Run full pipeline and compare with MATLAB outputs (1-to-1 matching)."""
        from tests.test_integration.matlab_compare import compare_isjc2_outputs, load_matlab_isjc2
        
        matlab_output_dir = isolate_srr25242877["matlab_output"]
        
        # Check MATLAB output exists (load_matlab_isjc2 handles REQUIRE_MATLAB_FILES)
        matlab_isjc2 = load_matlab_isjc2(matlab_output_dir)
        if matlab_isjc2 is None:
            pytest.skip("MATLAB reference output not available")
        
        # Create config (without ancestor for this test)
        config = self._create_config(isolate_srr25242877, cleared_output_dir)
        
        # Setup pipeline
        run_dir = self._setup_pipeline(config, return_run_dir=True)
        
        # Run full pipeline with step-by-step comparison
        pipeline = TestPipeline(config, matlab_output_dir)
        result = pipeline.run()
        
        # Check if pipeline produced any candidates
        assert len(result) > 0, \
            f"Pipeline found 0 candidates. Expected to match MATLAB output with {len(matlab_isjc2)} junctions. " \
            f"Check pipeline parameters and matching logic. Run dir: {run_dir}"
        
        isjc2_file = run_dir / "ISJC2.csv"
        assert isjc2_file.exists(), f"ISJC2.csv not found in {run_dir}"
        candidate_file = run_dir / "candidate_amplifications.csv"
        assert candidate_file.exists(), f"candidate_amplifications.csv not found in {run_dir}"
        
        # Load and compare outputs
        python_isjc2 = pd.read_csv(run_dir / "ISJC2.csv")
        print(f"\n=== Final Comparison ===")
        print(f"Python: {len(python_isjc2)} candidates")
        print(f"MATLAB: {len(matlab_isjc2)} candidates")
        
        # Compare outputs (1-to-1 matching required)
        compare_isjc2_outputs(python_isjc2, matlab_isjc2)
        print("✓ Comparison passed: 1-to-1 match with MATLAB")
    
    def test_isolate_pipeline_steps_without_ancestor(self, isolate_srr25242877, tmp_path, cleared_output_dir):
        """Test isolate pipeline step-by-step without ancestor (no MATLAB comparison)."""
        matlab_output_dir = isolate_srr25242877["matlab_output"]
        config = self._create_config(isolate_srr25242877, cleared_output_dir)
        
        # Setup pipeline
        self._setup_pipeline(config)
        
        print("\n=== Testing Isolate Pipeline (no ancestor) ===")
        pipeline = TestPipeline(config, matlab_output_dir)
        result = pipeline.run()
        
        print(f"\nFinal result: {len(result)} candidates")
        return result
    
    def test_ancestor_as_isolate_pipeline_steps(self, isolate_srr25242906, tmp_path, cleared_output_dir):
        """Test ancestor run as isolate (self-ancestor) - step-by-step comparison."""
        matlab_output_dir = isolate_srr25242906["matlab_output"]
        config = self._create_config(isolate_srr25242906, cleared_output_dir)
        config.anc_name = isolate_srr25242906["iso_name"]  # Self-ancestor
        
        # Setup pipeline
        self._setup_pipeline(config)
        
        print("\n=== Testing Ancestor Pipeline ===")
        pipeline = TestPipeline(config, matlab_output_dir)
        result = pipeline.run()
        
        print(f"\nFinal result: {len(result)} candidates")
        return result
    
    def test_isolate_pipeline_with_ancestor_steps(self, isolate_srr25242877, isolate_srr25242906, tmp_path, cleared_output_dir):
        """Test isolate pipeline with ancestor - step-by-step comparison with MATLAB."""
        from tests.test_integration.matlab_compare import compare_isjc2_outputs, load_matlab_isjc2
        
        matlab_output_dir = isolate_srr25242877["matlab_output"]
        
        # Check MATLAB output exists (load_matlab_isjc2 handles REQUIRE_MATLAB_FILES)
        matlab_isjc2 = load_matlab_isjc2(matlab_output_dir)
        if matlab_isjc2 is None:
            pytest.skip("MATLAB reference output not available")
        
        config = self._create_config(isolate_srr25242877, cleared_output_dir, anc_isolate=isolate_srr25242906)
        
        # Setup pipeline
        run_dir = self._setup_pipeline(config, return_run_dir=True)
        
        print("\n=== Testing Isolate Pipeline (with ancestor) ===")
        print("Note: Ancestor pipeline will run automatically if needed")
        pipeline = TestPipeline(config, matlab_output_dir)
        result = pipeline.run()  # Will run ancestor automatically if needed
        
        print(f"\nFinal result: {len(result)} candidates")
        
        # Verify output file exists
        python_isjc2_file = run_dir / "ISJC2.csv"
        assert python_isjc2_file.exists(), f"ISJC2.csv not found in {run_dir}"
        
        # Load and compare outputs
        python_isjc2 = pd.read_csv(python_isjc2_file)
        print(f"\n=== Comparing with MATLAB ISJC2.xlsx ===")
        print(f"Python: {len(python_isjc2)} candidates")
        print(f"MATLAB: {len(matlab_isjc2)} candidates")
        
        # Compare outputs (1-to-1 matching required)
        compare_isjc2_outputs(python_isjc2, matlab_isjc2)
        print("✓ Comparison passed: 1-to-1 match with MATLAB")
        
        return result
