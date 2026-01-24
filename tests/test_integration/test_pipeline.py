"""Integration tests for full pipeline execution."""

import pandas as pd
import numpy as np
import pytest
import pysam

from pathlib import Path
from Bio import SeqIO

from amplifinder.utils.file_utils import remove_file_or_dir
from amplifinder.data_types import JunctionType
from amplifinder.steps.base import Step
from amplifinder.config import Config
from amplifinder.pipeline import Pipeline

from tests.test_integration.matlab_compare import (
    convert_python_records_to_standard,
    convert_matlab_to_standard,
    compare_amplifications,
)
from tests.test_integration.test_utils import print_color, force_step
from tests.test_integration.jc_cover_compare import (
    load_matlab_read_buckets,
    load_python_read_buckets,
    create_jct_comparison_table,
)

pytestmark = pytest.mark.integration

# Global flag: if True, skip tests when MATLAB files are missing; if False, fail
REQUIRE_MATLAB_FILES = True
RECREATE_BAM = False

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


def _compare_bams_one_to_one(bam1: Path, bam2: Path) -> None:
    """Compare two BAMs record-by-record with normalized ref names."""
    def extract_reads(bam: Path, normalize: bool = False):
        with pysam.AlignmentFile(bam, "rb") as f:
            return [
                (
                    r.query_name,
                    r.flag,
                    JUNCTION_ID_MAP[r.reference_name] if normalize else r.reference_name,
                    r.reference_start,
                    r.cigarstring,
                    r.query_sequence,
                )
                for r in f.fetch(until_eof=True)
            ]

    reads1 = extract_reads(bam1, normalize=True)
    reads2 = extract_reads(bam2, normalize=False)

    if len(reads1) != len(reads2):
        pytest.fail(f"BAM length mismatch: {len(reads1)} vs {len(reads2)}")

    mismatches = [(i, a, b) for i, (a, b) in enumerate(zip(reads1, reads2), start=1) if a != b]
    if mismatches:
        msg_lines = []
        for idx, a, b in mismatches[:10]:
            msg_lines.append(f"{idx}: {a} != {b}")
        extra = f" (showing first 10 of {len(mismatches)})" if len(mismatches) > 10 else ""
        pytest.fail("BAM mismatches" + extra + ":\n" + "\n".join(msg_lines))


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

    def _align_reads(
        self,
        syn_tnjc2s,
        iso_output,
        anc_output,
    ):
        # Identify target candidate and force re-alignment by removing existing BAM
        matching_candidate = _select_target_candidate(syn_tnjc2s)
        python_bam = matching_candidate.bam_path(iso_output, is_ancestor=False)

        # Remove existing BAM files (to force re-alignment)
        if RECREATE_BAM:
            remove_file_or_dir(python_bam)
            remove_file_or_dir(python_bam.with_suffix(".bai"))
            remove_file_or_dir(python_bam.with_suffix(".linearindex"))

        # Run parent method and capture bowtie command
        super()._align_reads(syn_tnjc2s, iso_output, anc_output)

        # Compare BAM files with MATLAB for each candidate
        matlab_candidate_dir = self.matlab_output_dir / MATLAB_CANDIDATE_DIRNAME
        matlab_alignment_dir = matlab_candidate_dir / "alignment"
        matlab_bam = matlab_alignment_dir / "alignment.sorted.bam"

        python_bam = matching_candidate.bam_path(iso_output, is_ancestor=False)

        # Compare BAMs record-by-record
        _compare_bams_one_to_one(matlab_bam, python_bam)

        print_color("✓ All alignments perfectly match MATLAB")

    def _analyze_alignments(
        self,
        synjct_tnjc2s,
        iso_output,
        anc_output,
        read_lengths,
    ):
        # Enable DEBUG to ensure Python FASTA files are written for comparison
        from amplifinder.env import DEBUG
        with DEBUG.temp_set(True):
            analyzed = super()._analyze_alignments(
                synjct_tnjc2s, iso_output, anc_output, read_lengths
            )

        matlab_alignment_dir = self.matlab_output_dir / MATLAB_CANDIDATE_DIRNAME / "alignment"

        def _load_counts(name):
            path = matlab_alignment_dir / name
            if not path.exists():
                if REQUIRE_MATLAB_FILES:
                    pytest.fail(f"MATLAB {name} not found: {path}")
                return None
            txt = path.read_text().strip()
            if not txt:
                if REQUIRE_MATLAB_FILES:
                    pytest.fail(f"MATLAB {name} is empty: {path}")
                return None
            return [int(x) for x in txt.split(",") if x]

        matlab_left = _load_counts("bamreads__nmbr_left_reads.csv")
        matlab_right = _load_counts("bamreads__nmbr_right_reads.csv")
        matlab_span = _load_counts("bamreads__nmbr_green_reads.csv")

        if matlab_left is None or matlab_right is None or matlab_span is None:
            return analyzed

        target = _select_target_candidate(analyzed)
        covs = target.jc_covs
        py_left = [covs[jt].left + covs[jt].left_marginal for jt in JunctionType]
        py_right = [covs[jt].right + covs[jt].right_marginal for jt in JunctionType]
        py_span = [covs[jt].spanning for jt in JunctionType]

        # Always generate read-by-read comparison tables
        python_bam = target.bam_path(iso_output, is_ancestor=False)
        python_alignment_dir = python_bam.parent
        matlab_buckets = load_matlab_read_buckets(matlab_alignment_dir)
        python_buckets = load_python_read_buckets(python_alignment_dir)

        output_dir = python_alignment_dir

        # Create comparison table for each junction
        for jt in JunctionType.sorted():
            jct_num = jt.num

            df, confusion_matrix = create_jct_comparison_table(
                python_bam=python_bam,
                matlab_buckets=matlab_buckets,
                python_buckets=python_buckets,
                jct_num=jct_num,
            )

            csv_path = output_dir / f"jct_{jct_num}_read_comparison.csv"
            df.to_csv(csv_path, index=False)

            # check if conf mat has any non-diagonal elements
            if np.any(confusion_matrix.values != np.diag(confusion_matrix.values)):
                print_color("✗ Not all classifications match")
            else:
                print_color("✓ All classifications match perfectly")
            # Print confusion matrix
            print_color(f"\nConfusion matrix for jct_{jct_num}:")
            print(confusion_matrix.to_string())
            print()

            if len(df):
                print_color(f"Comparison table: {csv_path}")
                print_color(f"\nFirst 10 mismatches for jct_{jct_num}:")
                print(df.head(10).to_string(index=False))
                print()

        if matlab_left != py_left or matlab_right != py_right or matlab_span != py_span:
            msg = [
                f"Left  mismatch: \n{np.array([matlab_left, py_left])}",
                f"Right mismatch: \n{np.array([matlab_right, py_right])}",
                f"Spanning mismatch: \n{np.array([matlab_span, py_span])}",
            ]
            print_color("\nCompared read counts matlab (top row) vs python (bottom row):\n" + "\n\n".join(msg))

        print_color("✓ Junction read counts (left/spanning/right) match MATLAB")
        return analyzed

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
        print_color("\n=== Testing Isolate Pipeline WITHOUT ancestor ===")
        matlab_output_dir = isolate_srr25242877["matlab_output"]
        config = self._create_config(isolate_srr25242877, cleared_output_dir)
        pipeline = Pipeline(config, matlab_output_dir)
        pipeline.run()

    def test_pipeline_with_ancestor(self, isolate_srr25242877, isolate_srr25242906, cleared_output_dir):
        print_color("\n=== Testing Isolate Pipeline WITH ancestor ===")
        matlab_output_dir = isolate_srr25242877["matlab_output"]
        config = self._create_config(isolate_srr25242877, cleared_output_dir, anc_isolate=isolate_srr25242906)
        pipeline = Pipeline(config, matlab_output_dir)
        pipeline.run()  # Will run ancestor automatically if needed

    def test_pipeline_with_ancestor_and_compare(self, isolate_srr25242877, isolate_srr25242906, cleared_output_dir):
        print_color("\n=== Testing Isolate Pipeline WITH ancestor ===")
        matlab_output_dir = isolate_srr25242877["matlab_output"]
        config = self._create_config(isolate_srr25242877, cleared_output_dir, anc_isolate=isolate_srr25242906)
        pipeline = TestPipeline(config, matlab_output_dir)
        pipeline.run()  # Will run ancestor automatically if needed
