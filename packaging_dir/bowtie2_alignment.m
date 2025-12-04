function bowtie2_alignment(ISJC2,iso,jc_file,alignment_directory)

n = 6; %multithreading
bowtie_prms.score_min1 = 0;
bowtie_prms.score_min2 = -0.25;
bowtie_prms.mismatch_penalty = [5,5];
bowtie_prms.alignment_mode = '--end-to-end';
bowtie_prms.num_hits = 100;

outFile1 = 'alignment.sam'; 
outFile2 = 'alignment.sorted.bam';

for i = 1:height(ISJC2)
    isjc2 = ISJC2(i,:);    
    %isolate
    outDir = fullfile(isjc2.iso_outpath,isjc2.analysis_directory,alignment_directory);    
    output_exists = exist(strcat(outDir, filesep, outFile1),'file') | ...
        exist(strcat(outDir, filesep, outFile2),'file') ;
    if ~output_exists
        refPath = fullfile(isjc2.iso_outpath,isjc2.analysis_directory,jc_file);
        buildName = jc_file(1:end-6);
        fastqPath = iso.fastq_path(1); %isoalte
        run_bowtie_line(refPath,buildName,fastqPath,outDir,bowtie_prms)
    end
    %ancestor
    outDir = fullfile(isjc2.anc_outpath,isjc2.analysis_directory,alignment_directory);
    output_exists = exist(strcat(outDir, filesep, outFile1),'file') | ...
        exist(strcat(outDir, filesep, outFile2),'file') ;
    if ~output_exists
        refPath = fullfile(isjc2.anc_outpath,isjc2.analysis_directory,jc_file);
        buildName = jc_file(1:end-6);
        fastqPath = iso.fastq_path(2); %ancestor 
        run_bowtie_line(refPath,buildName,fastqPath,outDir,bowtie_prms)
    end
end
