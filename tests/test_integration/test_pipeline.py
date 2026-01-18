"""Integration tests for full pipeline execution."""

import pandas as pd
import numpy as np
import pytest
import pysam

from contextlib import contextmanager
from pathlib import Path
from Bio import SeqIO

from amplifinder.utils.file_utils import remove_file_or_dir
from amplifinder.data_types.enums import JunctionType
from amplifinder.steps.base import Step
from amplifinder.config import Config
from amplifinder.pipeline import Pipeline

from tests.test_integration.matlab_compare import (
    convert_python_records_to_standard,
    convert_matlab_to_standard,
    compare_amplifications,
    load_matlab_isjc2
)
from tests.test_integration.test_utils import print_color, force_step

pytestmark = pytest.mark.integration

# Global flag: if True, skip tests when MATLAB files are missing; if False, fail
REQUIRE_MATLAB_FILES = True

TARGET_LEFT = 873161
TARGET_RIGHT = 890745
MATLAB_CANDIDATE_DIRNAME = "chr_U00096_873161F_890745R_IS_41R"

# Mapping from junction numeric IDs to names
JUNCTION_ID_MAP = {str(j.num): j.name for j in JunctionType}


def _select_target_candidate(syn_tnjc2s):
    """Pick the 873161-890745 candidate; fallback to first if not found."""
    for c in syn_tnjc2s.to_records():
        if c.left == TARGET_LEFT and c.right == TARGET_RIGHT:
            return c
    pytest.fail(f"No matching candidate found for positions {TARGET_LEFT}-{TARGET_RIGHT}")


class TestPipeline(Pipeline):
    """Pipeline subclass for testing with step-by-step MATLAB comparisons."""
    __test__ = False  # avoid pytest collection warning for helper class

    def __init__(self, config, matlab_output_dir):
        super().__init__(config)
        self.matlab_output_dir = matlab_output_dir
        self.ref_tns = None  # Store reference TNs for IS name mapping

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
        result = super()._locate_tns_in_reference(genome)
        self.ref_tns = result  # Store for IS name mapping
        return result

    def _create_synthetic_junctions(self, filtered_tnjc2s, genome, ref_tns,
                                    iso_output, anc_output, read_lengths):
        with force_step():
            syn_tnjc2s = super()._create_synthetic_junctions(
                filtered_tnjc2s, genome, ref_tns, iso_output, anc_output, read_lengths,
            )

        matlab_fasta = self.matlab_output_dir / MATLAB_CANDIDATE_DIRNAME / "jc.fasta"
        if not matlab_fasta.exists():
            if REQUIRE_MATLAB_FILES:
                pytest.fail(f"MATLAB jc.fasta not found: {matlab_fasta}")
            return syn_tnjc2s

        target = _select_target_candidate(syn_tnjc2s)
        python_fasta = target.fasta_path(iso_output, is_ancestor=False)

        matlab_seqs = {JUNCTION_ID_MAP[rec.id]: str(rec.seq) for rec in SeqIO.parse(matlab_fasta, "fasta")}
        python_seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(python_fasta, "fasta")}

        diffs = [k for k in matlab_seqs if python_seqs[k] != matlab_seqs[k]]
        if diffs:
            pytest.fail(f"Junction FASTA mismatch for: {', '.join(diffs)}")
        
        if list(matlab_seqs.keys()) != list(python_seqs.keys()):
            pytest.fail(f"Junction FASTA key-order mismatch: {list(matlab_seqs.keys())} != {list(python_seqs.keys())}")

        print_color("✓ Junction FASTA perfectly matches MATLAB")
        return syn_tnjc2s

    def _create_ref_tn_junctions(self, tn_loc, genome, iso_output):
        result = super()._create_ref_tn_junctions(tn_loc, genome, iso_output)
        return result

    def _run_breseq(self, genome):
        result = super()._run_breseq(genome)
        return result

    def _create_tnjcs(self, breseq_jc, ref_tnjc, genome, iso_output):
        result = super()._create_tnjcs(breseq_jc, ref_tnjc, genome, iso_output)

        # Ensure single-locus TNJCs cover all MATLAB ISJC2 sides
        matlab_df = self._load_matlab_isjc2()
        if matlab_df is not None and len(matlab_df) > 0:
            matlab_positions = pd.concat([
                matlab_df['Positions_in_chromosome_1'],
                matlab_df['Positions_in_chromosome_2']
            ]).dropna().astype(int).to_numpy()

            tnjc_positions = result.df["pos2"].to_numpy()
            closest_diffs = np.array([int(np.min(np.abs(matlab_positions - pos))) for pos in tnjc_positions])
            closest_diffs_matlab = np.array([int(np.min(np.abs(tnjc_positions - pos))) for pos in matlab_positions])

            # Assert perfect matching
            assert np.all(closest_diffs == 0), \
                f"Python TNJC positions not matching MATLAB: " \
                f"{np.sum(closest_diffs > 0)} mismatches"
            assert np.all(closest_diffs_matlab == 0), \
                f"MATLAB positions not found in Python: " \
                f"{np.sum(closest_diffs_matlab > 0)} missing"

        return result

    def _create_tnjc2s(self, tnjc, genome, iso_output):
        result = super()._create_tnjc2s(tnjc, genome, iso_output)

        # Compare RawTnJc2 with MATLAB ISJC2.xlsx (positions, length, IS elements)
        matlab_df = self._load_matlab_isjc2()
        if matlab_df is not None:
            print_color("\nStep 6: Comparing RawTnJc2 with MATLAB ISJC2")

            # Convert both to standard format
            records = result.to_records()
            ref_tns_dict = self.ref_tns.to_dict()
            python_std = convert_python_records_to_standard(records, ref_tns_dict)
            matlab_std = convert_matlab_to_standard(matlab_df)

            compare_amplifications(python_std, matlab_std)
            print_color("✓ Comparison complete")
        else:
            assert not REQUIRE_MATLAB_FILES, "MATLAB file missing but REQUIRE_MATLAB_FILES=True"
            print_color(f"Python={len(result)} candidates (MATLAB file not found)")

        return result

    def _calc_amplicon_coverage(self, tnjc2, genome, iso_output):
        result = super()._calc_amplicon_coverage(tnjc2, genome, iso_output)
        return result

    def _export(self, analyzed, genome, iso_output):
        super()._export(analyzed, genome, iso_output)

        # Note: Final comparison done in test methods
        # Skipping intermediate export comparison here


@pytest.mark.slow
class TestPipelineStepByStep:
    """Pipeline integration tests with step-by-step MATLAB comparison."""

    @staticmethod
    def _setup_pipeline(config, return_run_dir=False):
        """Common setup: enable verbose reporting."""
        Step.set_global_verbose(True)
        if return_run_dir:
            return config.iso_run_dir
        return None

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

    def test_pipeline_without_ancestor(self, isolate_srr25242877, cleared_output_dir):
        """Pipeline on isolate only - compare with MATLAB."""
        matlab_output_dir = isolate_srr25242877["matlab_output"]

        # Check MATLAB output exists
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

        isjc2_file = run_dir / "tnjc2_exported.csv"
        assert isjc2_file.exists(), f"tnjc2_exported.csv not found in {run_dir}"
        candidate_file = run_dir / "candidate_amplifications.csv"
        assert candidate_file.exists(), f"candidate_amplifications.csv not found in {run_dir}"

        print_color("\n=== Final ISJC2 Comparison ===")
        print_color("✓ Comparison done in _create_tnjc2 step")

    def test_pipeline_with_ancestor(
            self,
            isolate_srr25242877,
            isolate_srr25242906,
            cleared_output_dir):
        """Pipeline on isolate + ancestor - compare with MATLAB."""
        matlab_output_dir = isolate_srr25242877["matlab_output"]

        # Check MATLAB output exists
        matlab_isjc2 = load_matlab_isjc2(matlab_output_dir)
        if matlab_isjc2 is None:
            pytest.skip("MATLAB reference output not available")

        config = self._create_config(isolate_srr25242877, cleared_output_dir, anc_isolate=isolate_srr25242906)

        # Setup pipeline
        run_dir = self._setup_pipeline(config, return_run_dir=True)

        print_color("\n=== Testing Isolate Pipeline (with ancestor) ===")
        pipeline = TestPipeline(config, matlab_output_dir)
        result = pipeline.run()  # Will run ancestor automatically if needed

        print_color(f"\nFinal result: {len(result)} candidates")

        # Verify output file exists
        python_isjc2_file = run_dir / "tnjc2_exported.csv"
        assert python_isjc2_file.exists(), f"tnjc2_exported.csv not found in {run_dir}"

        print_color("\n=== Final Comparison ===")
        print_color("✓ Comparison done in _create_tnjc2 step")
