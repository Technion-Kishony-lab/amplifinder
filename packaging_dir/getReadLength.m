function getReadLength(inpth,outpth)
if ~exist(fullfile(outpth,'read_length.mat'),'file')
    path_to_fastq = [char(inpth), filesep, 'data'] ;
    fn = dir([path_to_fastq filesep '*.fastq']) ;
    fn = fn(1).name ;
    fq = fastqread([path_to_fastq filesep fn],'Blockread',[1,10]) ;
    read_length = cellfun(@length,{fq.Sequence}) ;
    save(fullfile(outpth, 'read_length.mat'),'read_length')
end