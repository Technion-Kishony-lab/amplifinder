function breseq2mat(iso,ref,iso_inpath,iso_outpath)

global prms
BRESEQ_OUTPUT_SIZE_THRS = prms.BRESEQ_OUTPUT_SIZE_THRS ; 

breseq_output_file = fullfile(iso_inpath,'output','output.gd') ;
output_size_file = fullfile(iso_outpath,'output_size.mat') ;

% counting number of lines on output file.
if ~exist(output_size_file,'file')
    n = count_lines(char(breseq_output_file));
    if ~isnan(n)
        save(output_size_file,'n')
    end
else
    load(output_size_file, 'n')
end

coverage_error = false ;
size_error = n > BRESEQ_OUTPUT_SIZE_THRS;
if size_error
    fprintf('warning: isolate %s removed. output.gd is %d lines long.\n',iso{1},n)
end
iso_outfile = fullfile(iso_outpath, 'MUT.mat') ; 
[coverage_error,h] = parseBreseqCOV(iso,ref{1},iso_inpath,iso_outpath) ;
if coverage_error
    fprintf('warning: isolate %s removed. coverage files is only %d lines long.\n',iso{1},n)
end
parseBreseqMUT(breseq_output_file,iso_outfile) ;
get_breseq_ver(iso_inpath,iso_outpath) ;

getReadLength(iso_inpath,iso_outpath) ;
getMappedBases(iso_inpath,iso_outpath) ;

if size_error || coverage_error
    logger('low_quality_alignment','ERROR')
end