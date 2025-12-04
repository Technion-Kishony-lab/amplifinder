function run_samtools_line(outDir,min_qlen)

global prms
samtools_path = prms.SAMTOOLS_PATH;
outDir = char(outDir);
system([samtools_path filesep 'samtools view -bS --min-qlen ' min_qlen ' -@ 6 -o ' outDir filesep 'alignment.bam ' outDir filesep 'alignment.sam']);
system([samtools_path filesep 'samtools sort -@ 6 ' outDir filesep 'alignment.bam -o ' outDir filesep 'alignment.sorted.bam']);
system([samtools_path filesep 'samtools index ' outDir filesep 'alignment.sorted.bam']);
system(['rm ' outDir filesep 'alignment.sam']);
system(['rm ' outDir filesep 'alignment.bam']);
end
