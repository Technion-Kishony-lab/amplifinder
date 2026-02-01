"""Read length detection and junction length calculation step."""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from amplifinder.steps.base import OutputStep
from amplifinder.logger import logger, c
from amplifinder.utils.fasta import get_read_length_stats, count_total_bases


JUNCTION_LENGTH_TOLERANCE: float = 0.10
SAMPLE_PER_FILE: int = 1000


@dataclass
class ReadLengths:
    """Read length and junction arm length results."""
    read_len_iso: int
    read_len_anc: Optional[int]
    jc_arm_len_iso: int
    jc_arm_len_anc: Optional[int]


class ReadLenStep(OutputStep[ReadLengths]):
    """Determine read lengths from FASTQ and calculate junction lengths."""
    NAME = "Determine read lengths"

    def __init__(
        self, *,
        iso_fastq_path: Path,
        anc_fastq_path: Optional[Path] = None,
        iso_read_length: Optional[int] = None,
        anc_read_length: Optional[int] = None,
        min_num_bases: Optional[int] = None,
    ):
        self.iso_fastq_path = iso_fastq_path
        self.anc_fastq_path = anc_fastq_path
        self._iso_read_length = iso_read_length
        self._anc_read_length = anc_read_length
        self.min_num_bases = min_num_bases
        super().__init__()

    @staticmethod
    def _calc_read_length(fastq_dir: Path, provided_length: Optional[int], sample_type: str) -> int:
        """Calculate read length from FASTQ files."""
        if provided_length is not None:
            return provided_length

        stats = get_read_length_stats(fastq_dir, sample_per_file=SAMPLE_PER_FILE)
        logger.info(f"{sample_type:<9}: {c(str(stats), 'cyan')}")
        if not stats.is_uniform:
            logger.warning(f"{sample_type} read length is not uniform - may affect coverage accuracy")
        return stats.max_length

    def _determine_junction_arm_lengths(
        self, iso_read_len: int, anc_read_len: Optional[int]
    ) -> tuple[int, Optional[int]]:
        """Determine junction arm length (2x read length) for iso/anc with 10% tolerance."""
        iso_jc_arm = 2 * iso_read_len
        if anc_read_len is None:
            logger.info(f"Junction arm length: no ancestor; using 2*iso_read_len = {iso_jc_arm} for isolate junctions")
            return iso_jc_arm, None

        delta = abs(iso_read_len - anc_read_len)
        threshold = int(anc_read_len * JUNCTION_LENGTH_TOLERANCE)

        def _log_decision(prefix: str, using_len: int, func: callable) -> None:
            msg = f"{prefix} iso={iso_read_len}, anc={anc_read_len}, diff={delta}, " \
                f"{JUNCTION_LENGTH_TOLERANCE}% threshold (anc)={threshold}; using arm_length={using_len}"
            func(msg)
        anc_jc_arm = 2 * anc_read_len
        if delta <= threshold:
            iso_jc_arm = anc_jc_arm
            _log_decision("isolate-ancestor read lengths within threshold;", iso_jc_arm, logger.info)
        else:
            _log_decision("isolate-ancestor read lengths exceed threshold;", iso_jc_arm, logger.warning)
        return iso_jc_arm, anc_jc_arm

    def _calculate_output(self) -> ReadLengths:
        """Calculate read lengths and junction arm lengths."""
        read_len_iso = self._calc_read_length(self.iso_fastq_path, self._iso_read_length, "isolate")

        read_len_anc = None
        if self.anc_fastq_path:
            read_len_anc = self._calc_read_length(self.anc_fastq_path, self._anc_read_length, "ancestor")

        # Check minimum bases if specified
        if self.min_num_bases is not None:
            total_bases_iso = count_total_bases(self.iso_fastq_path, read_len_iso)
            if total_bases_iso < self.min_num_bases:
                logger.warning(
                    f"Low isolate sequencing depth: {total_bases_iso:,} bases "
                    f"(min recommended: {self.min_num_bases:,})"
                )

            if self.anc_fastq_path:
                total_bases_anc = count_total_bases(self.anc_fastq_path, read_len_anc)
                if total_bases_anc < self.min_num_bases:
                    logger.warning(
                        f"Low ancestor sequencing depth: {total_bases_anc:,} bases "
                        f"(min recommended: {self.min_num_bases:,})"
                    )

        jc_arm_len_iso, jc_arm_len_anc = self._determine_junction_arm_lengths(read_len_iso, read_len_anc)

        return ReadLengths(
            read_len_iso=read_len_iso,
            read_len_anc=read_len_anc,
            jc_arm_len_iso=jc_arm_len_iso,
            jc_arm_len_anc=jc_arm_len_anc,
        )
