function ISJC2 = bamAnalysis(ISJC2,alignment_directory,jc_file)

bamname = 'alignment.sorted.bam';
outputfile1 = 'biomap.mat';
outputfile2 = 'bamreads.mat';

for i = 1:height(ISJC2)
    fastafile = fullfile(ISJC2.iso_outpath(i),ISJC2.analysis_directory(i),jc_file);
    f = fastaread(fastafile);
    seq_length = length(f(1).Sequence);    
    %isolate
    load (fullfile(ISJC2.iso_outpath(i), 'fastq_read_length'),"read_length")
    path = fullfile(ISJC2.iso_outpath(i),ISJC2.analysis_directory(i),alignment_directory);
    if ~exist(fullfile(path,outputfile1),"file")
        parse_reads_from_bam(path,bamname,outputfile1);
    end
    load (fullfile(path,outputfile1),"RdStart","RdRef","RdLength","RdFlag","RdFull")
    if isempty(RdStart)
        green_reads = 1:7; left_reads = 1:7; right_reads = 1:7;
    else
        [green_reads,left_reads,right_reads]= count_reads(fullfile(path,outputfile2),seq_length,read_length,RdStart,RdRef,RdLength,RdFlag,RdFull);
    end
    ISJC2.iso_jc_cov_green(i) = {green_reads}; ISJC2.iso_jc_cov_left(i) = {left_reads}; ISJC2.iso_jc_cov_right(i) = {right_reads};
    %ancestor
    load (fullfile(ISJC2.anc_outpath(i), 'fastq_read_length'),"read_length")
    path = fullfile(ISJC2.anc_outpath(i),ISJC2.analysis_directory(i),alignment_directory);
    if ~exist(fullfile(path,outputfile1),"file")
        parse_reads_from_bam(path,bamname,outputfile1);
    end
    load (fullfile(path,outputfile1),"RdStart","RdRef","RdLength","RdFlag","RdFull")
    [green_reads,left_reads,right_reads]= count_reads(fullfile(path,outputfile2),seq_length,read_length,RdStart,RdRef,RdLength,RdFlag,RdFull);
    ISJC2.anc_jc_cov_green(i) = {green_reads}; ISJC2.anc_jc_cov_left(i) = {left_reads}; ISJC2.anc_jc_cov_right(i) = {right_reads};
end

end