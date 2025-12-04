function run_breseq(iso,ref_path,NCBI)
% assumes all references are in the same path. (all iso.ref_path are identical)

ref = string(strsplit(iso.ref{1},','));
if ~iso.isfinder(1)
    refPaths{1} = [char(iso.ref_path(1)) '/genbank'];
    refFiles{1} = strcat(ref, '.gb');
else
    refPaths{1} = [char(iso.ref_path(1)) '/fasta'];
    refFiles{1} = strcat(ref, '.fasta');
end

fastqPaths = iso.fastq_path ; 
outPaths = iso.breseq_path;  

for i = 1:height(iso)
    if exist(strcat(outPaths{i},filesep,'output',filesep,'output.gd'))
        continue
    else
        logger('cleaning breseq target directory','MESSAGE')
        system(['rm -r ', char(outPaths{i})]) %cleaning breseq target directory
        l = write_docker_line(refPaths{1},refFiles(1),fastqPaths{i},outPaths{i}) ;
        system(l)
    end
end

end

function l = write_docker_line(refPath,refFiles,fastqPath,outPath)
breseq_path = compile_breseq_path(refPath,fastqPath,outPath);
ref = char(strcat("-r /ref/", strjoin(refFiles{1}, " -r /ref/")));
fastqfiles = dir([char(fastqPath) filesep '*.fastq*']);
if isempty(fastqfiles)
    error(['Missing fastq files for ' fastqFolder])
end
fastqfiles = {fastqfiles.name};
fastqfiles = strcat('/fastq/', fastqfiles);
fastqfiles = strjoin(fastqfiles, ' ');
l = [breseq_path ' breseq -j 4 -o /out ' ref ' ' fastqfiles];
end

function breseq_path = compile_breseq_path(refPath,fastqPath,outPath)
ref_path = [char(refPath) ':/ref/'];
fastq_path = [char(fastqPath) ':/fastq/'];
outpath = [char(outPath) ':/out/'];
breseq_path = ['docker run --rm -u root -i -v ' char(ref_path) ' -v ' char(fastq_path) ...
    ' -v ' char(outpath) ' ummidock/breseq:0.32.1'] ;
% originally -it, now -i to allow multithreading
end