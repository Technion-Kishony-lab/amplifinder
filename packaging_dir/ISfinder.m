function ISfinder(ref,ref_path)
% Find IS elements by BLAST against IS-finder db.

global prms baseDir
ISDB_PATH = prms.ISDB_PATH ; 

outputName = [char(ref_path) filesep 'ISfinder' filesep ref '_IS_loc.mat'] ;
if ~exist(outputName,'file')
    [~,~] = mkdir ([char(ref_path),'ISfinder']) ; 
    F = fastaread(fullfile(baseDir,ISDB_PATH)) ; 
    len = cellfun(@length,{F.Sequence})' ; 
    headers = {F.Header}' ; 
    
    blast_output = [char(ref_path) '/ISfinder/' ref '_blast.txt'];
    blast_results = run_blast(fullfile(baseDir,ISDB_PATH), ...
    [char(ref_path) '/fasta/' ref '.fasta'],blast_output,1e-4);

    IS_loc = blast2ISloc(blast_results,headers,len);
    save(outputName, 'IS_loc')
end

