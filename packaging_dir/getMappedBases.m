function getMappedBases(inpth,outpth)
if ~exist(fullfile(outpth,'mapped_bases.mat'),'file')
    mapped_bases = [] ; 
    fn = strcat(inpth, filesep, '05_alignment_correction', filesep, 'summary.json') ;
    fin = fopen(fn,"r") ; 
    while ~feof(fin) 
        l = fgetl(fin) ; 
        tstf = strfind(l,'bases_mapped_to_reference') ; 
        if ~isempty(tstf)
            tstf = strfind(l,'total_bases_mapped_to_reference') ; 
            if isempty(tstf)
                tstf1 = strfind(l,':') ;
                tstf2 = strfind(l,',') ;
                mapped_bases = [mapped_bases ; str2num(l(tstf1+2:tstf2-1))] ;
            end
        end
    end
    save(fullfile(outpth, 'mapped_bases.mat'),'mapped_bases')
end