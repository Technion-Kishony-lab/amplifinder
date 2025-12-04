function name = efetch_genbank(raw_reference_name,ref_path)
% creates:
% genomeDB/genebank/XXXXX.mat: matlab genbank
% genomeDB/genebank/XXXXX.GB: NCBI genbank
% genomeDB/genebank/XXXXXref.mat: key genome properties

base_path = fullfile(ref_path, 'genbank');
key_filename = fullfile(base_path, [raw_reference_name 'ref.mat']);
baseURL = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
if ~exist(key_filename,'file')
    mkdir(base_path)
    
    % Get genebank (the matlab way):
    g = getgenbank(raw_reference_name) ;
    name = string(g.LocusName) ;
    length = str2num(g.LocusSequenceLength)  ; 
    circular = strcmp(g.LocusTopology,'circular') ; 
    ref_props = table(name,length,circular) ;
    save(key_filename, "ref_props")
    save(fullfile(base_path, name{1}), 'g')
    
    % get the raw genebank (whcih contains IS elements not in the matlab
    % genebank):
    gb_filename = fullfile(base_path,[name{1} '.gb']);
    websave(gb_filename, [baseURL 'efetch.fcgi?db=nuccore&id=' raw_reference_name '&rettype=gb&retmode=text']);
else
    load(key_filename, 'ref_props')
end

name = ref_props.name{1} ; 
