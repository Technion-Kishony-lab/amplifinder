"""Integration tests for full pipeline execution."""

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
from amplifinder.data_types import RecordTypedDF, Junction
from amplifinder.pipeline import Pipeline

pytestmark = pytest.mark.integration


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
        self.is_ancestor_run = False
    
    def _load_matlab_isjc2(self):
        """Load MATLAB ISJC2.xlsx if available."""
        import pandas as pd
        isjc2_xlsx = self.matlab_output_dir / "ISJC2.xlsx"
        if isjc2_xlsx.exists():
            return pd.read_excel(isjc2_xlsx)
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
        """Step 5: Match junctions to TN elements - compare with MATLAB."""
        result = super()._create_tnjc(breseq_jc, ref_tn_jc, ref_tn_end_seqs, genome, iso_output)
        
        # Compare with MATLAB (only for isolate runs with ancestor, not ancestor runs)
        if not self.is_ancestor_run:
            df = self._load_matlab_isjc2()
            if df is not None:
                # ISJC2 has junction pairs, so ISJC should be roughly 2x (each pair has 2 junctions)
                matlab_count = len(df) * 2  # Rough estimate
                print(f"Step 5: MATLAB ISJC2.xlsx has {len(df)} junction pairs (Python={len(result)} TN-associated junctions)")
                print(f"Step 5: MATLAB≈{matlab_count} IS-associated junctions (estimated), Python={len(result)} TN-associated junctions")
            else:
                print(f"Step 5: Python={len(result)} TN-associated junctions (MATLAB file not found)")
        else:
            print(f"Step 5: Python={len(result)} TN-associated junctions (ancestor run, no MATLAB comparison)")
        
        return result
    
    def _create_tnjc2(self, tnjc, genome, iso_output):
        """Step 6: Combine junction pairs - compare with MATLAB."""
        result = super()._create_tnjc2(tnjc, genome, iso_output)
        
        # Compare with MATLAB (only for isolate runs with ancestor, not ancestor runs)
        if not self.is_ancestor_run:
            df = self._load_matlab_isjc2()
            if df is not None:
                print(f"Step 6: MATLAB={len(df)} junction pairs, Python={len(result)} junction pairs")
            else:
                print(f"Step 6: Python={len(result)} junction pairs (MATLAB file not found)")
        else:
            print(f"Step 6: Python={len(result)} junction pairs (ancestor run, no MATLAB comparison)")
        
        return result


@pytest.mark.slow
class TestPipelineStepByStep:
    """Pipeline integration tests with step-by-step MATLAB comparison."""
    
    @staticmethod
    def _setup_pipeline(config, clear_output_dir, return_run_dir=False):
        """Common setup: clear output dir and enable verbose reporting."""
        run_dir = clear_output_dir(config)
        Step.set_global_verbose(True)
        return run_dir if return_run_dir else None
    
    @staticmethod
    def _get_test_output_root(matlab_output_dir):
        """Get test output root directory next to MATLAB outputs."""
        return matlab_output_dir.parent.parent.parent / "python_outputs"
    
    @staticmethod
    def _create_config(isolate, anc_isolate=None, anc_name=None, test_output_root=None):
        """Create Config with common defaults."""
        if test_output_root is None:
            matlab_output_dir = isolate["matlab_output"]
            test_output_root = TestPipelineStepByStep._get_test_output_root(matlab_output_dir)
        
        config_kwargs = {
            "iso_path": isolate["fastq_path"],
            "ref_name": "U00096",
            "iso_name": isolate["iso_name"],
            "output_dir": test_output_root / "output",
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
    
    def test_full_pipeline_matches_matlab(self, isolate_srr25242877, clear_output_dir):
        """Run full pipeline and compare with MATLAB outputs (1-to-1 matching)."""
        from tests.test_integration.matlab_compare import compare_isjc2_outputs, load_matlab_isjc2
        
        matlab_output_dir = isolate_srr25242877["matlab_output"]
        
        # Check MATLAB output exists
        matlab_isjc2 = load_matlab_isjc2(matlab_output_dir)
        if matlab_isjc2 is None:
            pytest.skip("MATLAB reference output not available")
        
        # Create config (without ancestor for this test)
        config = self._create_config(isolate_srr25242877)
        
        # Setup pipeline
        run_dir = self._setup_pipeline(config, clear_output_dir, return_run_dir=True)
        
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
    
    def test_isolate_pipeline_steps_without_ancestor(self, isolate_srr25242877, tmp_path, clear_output_dir):
        """Test isolate pipeline step-by-step without ancestor (no MATLAB comparison)."""
        matlab_output_dir = isolate_srr25242877["matlab_output"]
        config = self._create_config(isolate_srr25242877)
        
        # Setup pipeline
        self._setup_pipeline(config, clear_output_dir)
        
        print("\n=== Testing Isolate Pipeline (no ancestor) ===")
        pipeline = TestPipeline(config, matlab_output_dir)
        result = pipeline.run()
        
        print(f"\nFinal result: {len(result)} candidates")
        return result
    
    def test_ancestor_as_isolate_pipeline_steps(self, isolate_srr25242906, tmp_path, clear_output_dir):
        """Test ancestor run as isolate (self-ancestor) - step-by-step comparison."""
        matlab_output_dir = isolate_srr25242906["matlab_output"]
        config = self._create_config(isolate_srr25242906)
        config.anc_name = isolate_srr25242906["iso_name"]  # Self-ancestor
        
        # Setup pipeline
        self._setup_pipeline(config, clear_output_dir)
        
        print("\n=== Testing Ancestor Pipeline ===")
        pipeline = TestPipeline(config, matlab_output_dir)
        pipeline.is_ancestor_run = True
        result = pipeline.run()
        
        print(f"\nFinal result: {len(result)} candidates")
        return result
    
    def test_isolate_pipeline_with_ancestor_steps(self, isolate_srr25242877, isolate_srr25242906, tmp_path, clear_output_dir):
        """Test isolate pipeline with ancestor - step-by-step comparison with MATLAB."""
        from tests.test_integration.matlab_compare import compare_isjc2_outputs, load_matlab_isjc2
        
        matlab_output_dir = isolate_srr25242877["matlab_output"]
        
        # Check MATLAB output exists
        matlab_isjc2 = load_matlab_isjc2(matlab_output_dir)
        if matlab_isjc2 is None:
            pytest.skip("MATLAB reference output not available")
        
        config = self._create_config(isolate_srr25242877, anc_isolate=isolate_srr25242906)
        
        # Setup pipeline
        run_dir = self._setup_pipeline(config, clear_output_dir, return_run_dir=True)
        
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
