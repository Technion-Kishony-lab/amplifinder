function samtools_analysis(ISJC2,alignment_directory)

outFile = 'alignment.sorted.bam'; 

for i = 1:height(ISJC2)
    isjc2 = ISJC2(i,:);    
    %isolate
    outDir = fullfile(isjc2.iso_outpath,isjc2.analysis_directory,alignment_directory);    
    output_exists = exist(strcat(outDir, filesep, outFile),'file');
    if ~output_exists
        load (fullfile(isjc2.iso_outpath, 'fastq_read_length'),"read_length")
        min_qlen = num2str(round(read_length*0.9));
        run_samtools_line(outDir,min_qlen)
    end
    %ancestor
    outDir = fullfile(isjc2.anc_outpath,isjc2.analysis_directory,alignment_directory);
    output_exists = exist(strcat(outDir, filesep, outFile),'file');
    if ~output_exists
        load (fullfile(isjc2.anc_outpath, 'fastq_read_length'),"read_length")
        min_qlen = num2str(round(read_length*0.9));
        run_samtools_line(outDir,min_qlen)
    end
end
