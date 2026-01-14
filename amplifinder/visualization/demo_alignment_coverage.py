"""Demo script for alignment coverage plotting with dummy data."""
import random

from amplifinder.data_types import JunctionReadCounts, JunctionType
from amplifinder.utils.timing import timer
from amplifinder.visualization.plot_alignments import plot_junctions_coverage


def main():
    random.seed(42)
    
    jct_lengths = {jt: 600 for jt in JunctionType}
    read_len = 100
    min_overlap_len = 12
    
    # Generate dummy reads for each junction type
    alignment_data = {}
    alignment_data_anc = {}
    jc_covs = {}
    jc_covs_anc = {}
    jc_calls = {}
    jc_calls_anc = {}
    
    for jt in JunctionType:
        length = jct_lengths[jt]
        arm_len = length // 2
        
        # Generate dummy reads
        jc_cov_iso = JunctionReadCounts()
        jc_cov_anc = JunctionReadCounts()
        reads_iso = []
        reads_anc = []
        
        # Isolate: 255 reads
        for _ in range(255):
            start = random.randint(0, length - read_len)
            end = min(start + read_len + random.randint(-10, 10), length)
            read_type = jc_cov_iso.add_read(start, end, arm_len, min_overlap_len)
            reads_iso.append((start, end, read_type))
        
        # Ancestor: 130 reads
        for _ in range(130):
            start = random.randint(0, length - read_len)
            end = min(start + read_len + random.randint(-10, 10), length)
            read_type = jc_cov_anc.add_read(start, end, arm_len, min_overlap_len)
            reads_anc.append((start, end, read_type))
        
        alignment_data[jt] = reads_iso
        jc_covs[jt] = jc_cov_iso
        jc_calls[jt] = random.choice([True, False, None])
        
        alignment_data_anc[jt] = reads_anc
        jc_covs_anc[jt] = jc_cov_anc
        jc_calls_anc[jt] = random.choice([True, False, None])
    
    with timer('plot_junctions_coverage'):
        plot_junctions_coverage(
            alignment_data=alignment_data,
            alignment_data_anc=alignment_data_anc,
            jc_lengths=jct_lengths,
            jc_covs=jc_covs,
            jc_covs_anc=jc_covs_anc,
            jc_calls=jc_calls,
            jc_calls_anc=jc_calls_anc,
            title='Alignment Coverage Demo (Dummy Reads)',
            output_path='demo_alignment_coverage.png',
            min_overlap_len=min_overlap_len,
            read_len=read_len,
        )


if __name__ == '__main__':
    main()
