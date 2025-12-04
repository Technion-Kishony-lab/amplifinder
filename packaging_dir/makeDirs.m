function makeDirs(iso)
ref_paths = cellfun(@(x)split(x,','), iso.ref_path, 'UniformOutput',false);
u = unique(cat(1,ref_paths{:}));
for i = 1:length(u)
    [~,~] = mkdir(char(u(i)));
    [~,~] = mkdir([char(u(i)) '/fasta']);
    [~,~] = mkdir([char(u(i)) '/genbank']);
    [~,~] = mkdir([char(u(i)) '/ISfinder']);
end
[~,~] = mkdir('BRESEQ');
[~,~] = mkdir('figures');
[~,~] = mkdir('output');

%make isolate specific dirs
mkdir(char(iso.iso_outpath(1))); mkdir(char(iso.iso_outpath(2))); 