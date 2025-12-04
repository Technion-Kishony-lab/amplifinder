function AmpliFinder(varargin)
% inputs: 
% -iso_path: path to isolate [fastq / fastq.gz] files.
% -iso_breseq_path: path to isolate breseq output (default = "BRESEQ"+name of isolate).
% -iso_name: name of isolate (used to name output directories, default="isolate"). 
% -anc_path: path to ancestral isolate [fastq / fastq.gz] files. 
% -anc_name = name of ancestor (used to name output directories, default="ancestor").
% -anc_breseq_path: path to ancestor breseq output (default = "BRESEQ"+name of isolate).
% -ref_path: path to reference genome files. 
% -ref_name: name of reference genome. (for '-ref_name CP0000', a 'fasta/CP0000.fasta' file is expected). 
% -steps: either "pre" or "all". "pre" referes to ...(default="all")
% -NCBI: use ncbi reference with ref_name (default="true")
% -ISfinder: identify reference IS elements based on ISfinder. a 'true' value is required for non-NCBI genomes (default="false")
% a 'config.txt' file is expected in 'installation directory' and can be overridden by a file with the same name in the pwd. 

warning('OFF', 'MATLAB:table:PreallocateCharWarning')
warning('OFF', 'MATLAB:table:RowsAddedExistingVars')
warning('OFF', 'MATLAB:MKDIR:DirectoryExists');
warning('OFF', 'MATLAB:DELETE:FileNotFound');

global prms baseDir workingDir logDir

workingDir = pwd;

if isdeployed
    % For deployed application, use ctfroot
    baseDir = [ctfroot '/AmpliFinder'];
else
    % For MATLAB environment, use the script directory
    baseDir = fileparts(mfilename('fullpath'));
end

% Add subfunctions directory
if ~isdeployed
    subfunctionsDir = fullfile(baseDir);
    addpath(genpath(subfunctionsDir));
end

[iso,NCBI,steps] = parse_AmpliFinder_arguments(varargin);
prms = readConfig;
makeDirs(iso) %making required directories
logDir = char(iso.iso_outpath(1)); 
logger(sprintf('Starting analysis of %s',iso.iso(1)),"MESSAGE")
%% take read length stats and uniformity 
iso = analyze_fastq_files(iso); 

%% confirm isolate and ancestor sequencing compatibility
if any(iso.read_length./iso.read_length'>2)
    logger('Read length of isolate and ancestor is incompatible', 'ERROR')
end

%% Get reference genome from ncbi / local
% by default, create genomeDB folder with all relevant genome files (see get_reference):
iso = curate_reference(iso,NCBI);

%% BRESEQ
run_breseq(iso,NCBI)

%% parse breseq output into .mat
breseq2mat(iso.iso(1),iso.ref(1),iso.breseq_path(1),iso.iso_outpath{1}) ;   %isolate
breseq2mat(iso.iso(2),iso.ref(2),iso.breseq_path(2),iso.iso_outpath{2}) ;   %ancestor

%% analyze pairs of junctions
%loads JC, COV, IS_loc
xlsx_output_filename = [iso.iso_outpath{1} filesep 'ISJC2.xlsx'];

if ~exist(xlsx_output_filename,"file")
    ISJC2 = curate_candidate_amplicons(iso(1,:),iso.iso_outpath(2));
    export_ISJC2(ISJC2,iso);
end
load(fullfile(iso.iso_outpath(1),'ISJC2.mat'),"ISJC2")

%% Write synthetic junctions to fasta

% identify transpositions and classify ISJC2 accordingly 
ISJC2 = classify_ISJC2(ISJC2);

%filter for high coverage candidate amplicons OR deletions
z = (ISJC2.amplicon_coverage>prms.COPY_NMBR_THRS | ISJC2.amplicon_coverage<prms.DEL_COPY_NMBR_THRS) & ISJC2.amplicon_length>prms.FILTER_AMPLICON_LENGTH;
filteredISJC2 = ISJC2(z,:);

if height(filteredISJC2)>0 && strcmp(steps,"noAlignment")
    logger(sprintf('Completed partial analysis of %s. Skipped alignment to synthetic junctions.',iso.iso(1)),"MESSAGE")
    return
end

%writing fasta of filtered
plus_anc = true;
jc_file = 'jc.fasta';
junction2fasta(filteredISJC2,plus_anc,jc_file)

%% Running bowtie2 + samtools on synthetic junctions 

alignment_directory = 'alignment';
for i = 1:height(filteredISJC2)
    mkdir (fullfile(filteredISJC2.iso_outpath(i),filteredISJC2.analysis_directory(i),alignment_directory))
    mkdir (fullfile(filteredISJC2.anc_outpath(i),filteredISJC2.analysis_directory(i),alignment_directory))
end

bowtie2_alignment(filteredISJC2,iso,jc_file,alignment_directory)
samtools_analysis(filteredISJC2,alignment_directory)

%% BAM analysis

filteredISJC2 = bamAnalysis(filteredISJC2,alignment_directory,jc_file);
filteredISJC2 = classify_candidates(filteredISJC2);

%% export of final analysis output

export_filteredISJC2(filteredISJC2,iso);
logger(sprintf('Completed analysis of %s',iso.iso(1)),"MESSAGE")
