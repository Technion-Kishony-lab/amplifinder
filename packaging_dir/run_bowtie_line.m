function run_bowtie_line(refPath,buildName,fastqPath,outDir,bowtie_prms)

score_min1 = bowtie_prms.score_min1;
score_min2 = bowtie_prms.score_min2;
mismatch_penalty = bowtie_prms.mismatch_penalty;
alignment_mode = bowtie_prms.alignment_mode;
num_hits = bowtie_prms.num_hits;

ref = char(refPath);
outDir = char(outDir);
fastqfiles = dir([char(fastqPath) filesep '*.fastq*']);
if isempty(fastqfiles)
    error(['Missing fastq files for ' fastqPath])
end
fastqfiles = {fastqfiles.name};
fastqfiles = strcat(repmat(strcat(fastqPath, filesep),size(fastqfiles)), fastqfiles); 
fastqfiles = strjoin(fastqfiles, ',');

mismatch_penalty_text = sprintf('--mp %d,%d',mismatch_penalty(1),mismatch_penalty(2));
if strcmp(alignment_mode,'--local')
    score_min_text = sprintf('--score-min G,%d,%d', score_min1, score_min2);
else
    score_min_text = sprintf('--score-min L,%d,%1.1f', score_min1, score_min2);
end
num_hits_text = num2str(num_hits);

system(['bowtie2-build -f ' ref ' ' outDir filesep buildName]);
system(['bowtie2 -p 4 ' alignment_mode ' -k ' num_hits_text ' --reorder ' mismatch_penalty_text ' ' score_min_text ' -x ' outDir filesep buildName ' -U ' char(fastqfiles) ' -S ' outDir filesep 'alignment.sam']) ;

end
