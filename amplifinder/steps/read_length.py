"""Read length detection and junction length calculation step."""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from amplifinder.steps.base import OutputStep
from amplifinder.logger import info, warning
from amplifinder.utils.fasta import get_read_length_stats


JUNCTION_LENGTH_TOLERANCE: float = 0.10
SAMPLE_PER_FILE: int = 1000


@dataclass
class ReadLengths:
    """Read length and junction length results."""
    iso_read_length: int
    anc_read_length: Optional[int]
    iso_junction_length: int
    anc_junction_length: Optional[int]


class ReadLenStep(OutputStep[ReadLengths]):
    """Determine read lengths from FASTQ and calculate junction lengths."""

    def __init__(
        self, *,
        iso_fastq_path: Path,
        anc_fastq_path: Optional[Path] = None,
        iso_read_length: Optional[int] = None,
        anc_read_length: Optional[int] = None,
    ):
        self.iso_fastq_path = iso_fastq_path
        self.anc_fastq_path = anc_fastq_path
        self._iso_read_length = iso_read_length
        self._anc_read_length = anc_read_length
        super().__init__()

    @staticmethod
    def _calc_read_length(fastq_dir: Path, provided_length: Optional[int], sample_type: str) -> int:
        """Calculate read length from FASTQ files."""
        if provided_length is not None:
            return provided_length
        
        stats = get_read_length_stats(fastq_dir, sample_per_file=SAMPLE_PER_FILE)
        info(f"{sample_type:<9}: {stats}")
        if not stats.is_uniform:
            warning(f"{sample_type} read length is not uniform - may affect junction coverage accuracy")
        return stats.max_length

    def _determine_junction_lengths(self, iso_len: int, anc_len: Optional[int]) -> tuple[int, Optional[int]]:
        """Determine junction flank length for iso/anc with 10% tolerance."""
        if anc_len is None:
            info(f"Junction length: no ancestor; using isolate length {iso_len} for isolate junctions")
            return iso_len, None

        delta = abs(iso_len - anc_len)
        threshold = int(anc_len * JUNCTION_LENGTH_TOLERANCE)

        def _log_decision(prefix: str, using_len: int, func: callable) -> None:
            msg = f"{prefix} iso={iso_len}, anc={anc_len}, diff={delta}, " \
                f"{JUNCTION_LENGTH_TOLERANCE}% threshold (anc)={threshold}; using {using_len}"
            func(msg)

        if delta <= threshold:
            iso_jc_len = anc_len
            _log_decision("isolate-ancestor read lengths within threshold;", iso_jc_len, info)
        else:
            iso_jc_len = iso_len
            _log_decision("isolate-ancestor read lengths exceed threshold;", iso_jc_len, warning)    
        return iso_jc_len, anc_len

    def _calculate_output(self) -> ReadLengths:
        """Calculate read lengths and junction lengths."""
        iso_read_len = self._calc_read_length(self.iso_fastq_path, self._iso_read_length, "isolate")
        
        anc_read_len = None
        if self.anc_fastq_path:
            anc_read_len = self._calc_read_length(self.anc_fastq_path, self._anc_read_length, "ancestor")
        
        iso_jc_len, anc_jc_len = self._determine_junction_lengths(iso_read_len, anc_read_len)
        
        return ReadLengths(
            iso_read_length=iso_read_len,
            anc_read_length=anc_read_len,
            iso_junction_length=iso_jc_len,
            anc_junction_length=anc_jc_len,
        )
