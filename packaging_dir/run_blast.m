function t = run_blast(database,query,output_file,evalue,num_alignments,output_format)
% run blast.
% returns a table with blast results.

global prms baseDir
BLASTN_PATH = prms.BLASTN_PATH;

if nargin<4
    evalue = 1e-4;
end
if nargin<6
    output_format = '10';
end

if nargin<5
    num_alignments = 10000;
end
if ~exist(output_file,'file')
    if isequal(output_format,'-dust no')
        eval(['!' BLASTN_PATH 'blastn -db ' database ' -query ' query ' -evalue ' num2str(evalue) ...
        ' -num_alignments ' num2str(num_alignments) ' ' output_format ' > ' output_file])
    else
        eval(['!' BLASTN_PATH 'blastn -db ' database ' -query ' query ' -evalue ' num2str(evalue) ...
        ' -num_alignments ' num2str(num_alignments) ' -outfmt ' output_format ' > ' output_file])
    % default for 10: qaccver saccver pident length mismatch gapopen qtsart qend sstart send evalue bitscore'
    end
end
fprintf('loading opts from %s\n',[baseDir '/opts'])
load ([baseDir '/opts'], "opts")
if isequal(output_format,'10')
    t = readtable(output_file,opts) ;
    t.Properties.VariableNames = {'query','subject','percent_identical','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore'} ;
else
    t = [];
end
end
