"""Demo script for alignment coverage plotting with dummy data."""
import random

from amplifinder.config import AlignmentClassifyParams
from amplifinder.data_types import JunctionReadCounts, JunctionType
from amplifinder.steps.jct_coverage.alignment_data import SingleAlignment
from amplifinder.steps.jct_coverage.cigar import Cigar, merge_consecutive_cigar_ops
from amplifinder.steps.jct_coverage.read_type import get_hit_type
from amplifinder.utils.timing import timer
from amplifinder.visualization.plot_jc_alignments import plot_jc_alignments


def generate_dummy_hits(n_reads, length, read_len, arm_len, align_params, indel_at_junction=None):
    """Generate dummy alignment reads.

    Args:
        indel_at_junction: None, 'deletion', or 'insertion' to add indel at junction
    """
    jc_cov = JunctionReadCounts()
    reads = []
    for i in range(n_reads):
        start = random.randint(0, length - read_len)
        total_len = read_len + random.randint(-10, 10)

        # Generate cigar with SNPs
        cigar = Cigar()
        for pos in range(total_len):
            jc_pos = start + pos - arm_len
            if jc_pos != 0:
                prob_snp = 0.95 if jc_pos == -50 else 0.05
                is_snp = random.random() < prob_snp
                cigar.append((8, 1) if is_snp else (7, 1))
            elif indel_at_junction == 'deletion':
                cigar.append((2, 5))
            elif indel_at_junction == 'insertion':
                cigar.append((1, 5))
            else:
                cigar.append((7, 1))

        cigar = merge_consecutive_cigar_ops(cigar)
        end = start + cigar.get_total_length()
        read_type = get_hit_type(start + 1, end, arm_len, align_params)
        if read_type is not None:
            alignment = SingleAlignment(
                read_type=read_type, start=start, end=end, cigar=cigar, bam_index=i
            )
            reads.append(alignment)
            jc_cov.increment(read_type)

    return reads, jc_cov


def main():
    random.seed(42)

    jc_arm_len = 300
    read_len = 100
    min_overlap_len = 12
    align_params = AlignmentClassifyParams(min_overlap_len=min_overlap_len)

    # Generate dummy reads for each junction type
    alignment_data = {}
    alignment_data_anc = {}
    jc_covs = {}
    jc_covs_anc = {}
    jc_calls = {}
    jc_calls_anc = {}

    for jt in JunctionType:
        jct_length = jc_arm_len * 2

        hits_iso, jc_cov_iso = generate_dummy_hits(
            255, jct_length, read_len, jc_arm_len, align_params, indel_at_junction='insertion'
        )
        hits_anc, jc_cov_anc = generate_dummy_hits(
            130, jct_length, read_len, jc_arm_len, align_params, indel_at_junction='deletion'
        )

        alignment_data[jt] = hits_iso
        jc_covs[jt] = jc_cov_iso
        jc_calls[jt] = random.choice([True, False, None])

        alignment_data_anc[jt] = hits_anc
        jc_covs_anc[jt] = jc_cov_anc
        jc_calls_anc[jt] = random.choice([True, False, None])

    with timer('plot_junctions_coverage'):
        plot_jc_alignments(
            alignment_data=alignment_data,
            alignment_data_anc=alignment_data_anc,
            jc_covs=jc_covs,
            jc_covs_anc=jc_covs_anc,
            jc_calls=jc_calls,
            jc_calls_anc=jc_calls_anc,
            jc_arm_len_iso=jc_arm_len,
            jc_arm_len_anc=jc_arm_len,
            read_len_iso=read_len,
            read_len_anc=read_len,
            title='Alignment Coverage Demo (Dummy Reads)',
            output_path='demo_alignment_coverage.png',
            alignment_classify_params=AlignmentClassifyParams(min_overlap_len=min_overlap_len),
        )


if __name__ == '__main__':
    main()
