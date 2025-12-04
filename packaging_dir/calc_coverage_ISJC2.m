function [amplicon_copy_number_dist,ISJC2] = calc_coverage_ISJC2(ISJC2, ref_props, iso_outpath, anc_outpath)

global prms

h = height(ISJC2);
iso_coverage = nan(h,1);
anc_coverage = nan(h,1);
cov = cell(h,1);
anc_cov = cell(h,1);
median_genome_cov = zeros(h,1);
anc_median_genome_cov = zeros(h,1);
f = find(ISJC2.amplicon_length>prms.MIN_AMPLICON_LENGTH ...
    & ISJC2.amplicon_length<prms.MAX_AMPLICON_LENGTH)';

% seperate loops for isolate and ancestor to prevent reloading of coverage files in get_coverage_in_range
for i = f
    cov_file = strjoin([iso_outpath filesep ref_props.name(ISJC2.scaf_Chr(i)) '_COV.mat'], '');
    [cov{i}, median_genome_cov(i)] = get_coverage_in_range(...
        cov_file, ISJC2.pos_Chr(i,1), ISJC2.pos_Chr(i,2), ISJC2.span_origin(i));
    iso_coverage(i) = double(calc_average_coverage(cov{i})) ./ double(median_genome_cov(i));
end

for i = f
    cov_file = strjoin([anc_outpath filesep ref_props.name(ISJC2.scaf_Chr(i)) '_COV.mat'], '');
    [anc_cov{i}, anc_median_genome_cov(i)] = get_coverage_in_range(...
        cov_file, ISJC2.pos_Chr(i,1), ISJC2.pos_Chr(i,2), ISJC2.span_origin(i));
    anc_coverage(i) = double(calc_average_coverage(anc_cov{i})) ./ double(anc_median_genome_cov(i));
end

ISJC2.amplicon_coverage = iso_coverage ./ anc_coverage;

X1 = prms.NCP_LIMIT1;
X2 = prms.NCP_LIMIT2;
n = prms.NCP_N;
amplicon_copy_number_dist = zeros(height(ISJC2),n-1);
amplicon_coverage_mode = nan(height(ISJC2),1);
for i = f
    cp = double(cov{i}) ./ double(median_genome_cov(i));
    anc_cp = double(anc_cov{i}) ./ double(anc_median_genome_cov(i));
    ncp = cp ./(double(anc_cov{i})/double(anc_median_genome_cov(i)));
    ncp(ncp==Inf) = nan;
    
    ncp(ncp<10^X1) = 10^X1;
    ncp(ncp>10^X2) = 10^X2;
    hc_edges = logspace(X1,X2,n);
    hc = histcounts(ncp,logspace(X1,X2,n));
    amplicon_copy_number_dist(i,:) = hc;
    
    mx = max(hc(2:end));  % skipping first bin to avoid including deletions or coverage below 10^X1
    ff = find(hc==mx,1,'last');
    amplicon_coverage_mode(i) = mean(hc_edges(ff:ff+1));
end

ISJC2.amplicon_coverage_mode = amplicon_coverage_mode; 

end