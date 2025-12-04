close all
clear all
clc
warning('OFF', 'MATLAB:table:PreallocateCharWarning')
warning('OFF', 'MATLAB:table:RowsAddedExistingVars')
RUN_TIME = false;
global prms

%% general paths and parameters
addpath(genpath('scripts'))
[~,~] = mkdir('AmpliFinderWorkspace');

%% isolates
isoList = readIsolatesTable('isolates.xlsx'); 
fprintf('after filtering: analysing %d isolates\n\n', height(isoList)-1)

%% call AmpliFinder 
steps = 'all'; time_out_limit = 100;
wd = fullfile(pwd,'AmpliFinderWorkspace');
cd (wd)
ref_path = '/zdata/user-data/idan/small_projects/AmpliFinder/gDB1';

ampliFinder_path = '/zdata/user-data/idan/small_projects/AmpliFinder/v1.0/AmpliFinder1.0.38/for_redistribution_files_only'; 
runtime_path = '/zdata/user-data/idan/small_projects/AmpliFinder/v1.0/AmpliFinder1.0.29/installed/v911';

for i = 1:height(isoList)-1

        iso_name  = char(isoList.iso(i));
        iso_path = char(fullfile(isoList.iso_fastq_path(i),iso_name));
        iso_breseq_path = char(fullfile(isoList.breseq_output_path(i),iso_name));
        anc_name  = char(isoList.iso(end));
        anc_path = char(fullfile(isoList.iso_fastq_path(end),anc_name));
        anc_breseq_path = char(fullfile(isoList.breseq_output_path(end),anc_name));
        ref_name = char(isoList.input_ref(i));
        isfinder = char(string(isoList.isfinder(i)));
    
        fprintf('running AmpliFinder for isolate %s\n',iso_name)
        if RUN_TIME
            system(['sh ' ampliFinder_path filesep 'run_AmpliFinder.sh ' runtime_path...
                ' -iso_path ' iso_path ' -iso_breseq_path ' iso_breseq_path ...
                ' -iso_name ' iso_name ' -anc_path ' anc_path ...
                ' -anc_breseq_path ' anc_breseq_path ' -anc_name ' anc_name ...
                ' -ref_name ' ref_name ' -ref_path ' ref_path ' -ISfinder ' isfinder ' -steps ' steps]);
        else
            addpath(genpath('/zdata/user-data/rkishony/amplifinder/packaging_dir'))
            AmpliFinder ('-iso_path', iso_path, '-iso_breseq_path', iso_breseq_path, ...
                '-iso_name', iso_name, '-anc_path', anc_path, ...
                '-anc_breseq_path', anc_breseq_path, '-anc_name', anc_name, ...
                '-ref_name', ref_name, '-ref_path', ref_path, '-ISfinder', isfinder);
        end
end
cd ..

fprintf('AmpliFinder completed for %d isolates.\n',height(isoList)-1)

%% overview:
% *AmpliFinder outputs*
% each isolate has the following files as an AmpliFinder output:
% (1) classified_amplifications.mat
% (2) ISJC2.xlsx
% (3) XXXXX_COV.mat
% (4) ref_props.mat
% (5) IS_loc.mat
% (6) read_length.mat
%    Note: (1) is a filtered version of (2), which also include columns
%    summarizing the 2nd pass alignment analysis and a few changes in column
%    names. Filtering is by: copy number (median) of the amplicon and length.
%    Copy number is (>prms.COPY_NMBR_THRS) | (<prms.DEL_COPY_NMBR_THRS). 
%    Length is explicitly (>prms.FILTER_AMPLICON_LENGTH)
%    and implicitly (<prms.MAX_AMPLICON_LENGTH) as no coverage is determined for them. 
% For each entry/line in (1) a directory is generated in the both the isolate and the ancestor 
% directories. These directories contain the following AmpliFinder output:
% (a) bamreads.mat



