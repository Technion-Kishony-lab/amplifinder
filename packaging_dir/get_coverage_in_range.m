function [coverage, median_chr] = get_coverage_in_range(coverage_file, pos1, pos2, span_origin)
% returns the coverage in the specified genome between two positions

persistent previous_coverage_file cov previous_median_chr

if ~strcmp(previous_coverage_file, coverage_file)
    ld = load(coverage_file, 'COV');
    cov = ld.COV;
    median_chr = calc_average_coverage(cov);
    previous_median_chr = median_chr;
    previous_coverage_file = coverage_file;
else
    median_chr = previous_median_chr;
end

if span_origin
    coverage = cov([pos2:end, 1:pos1]);
else
    coverage = cov(pos1:pos2);
end
end