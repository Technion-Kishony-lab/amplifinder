function get_breseq_ver(inpth,outpth)

filename = fullfile(outpth,'breseq_ver.mat') ;
if ~exist(fullfile(outpth,'breseq_ver.mat'),'file')
    fid = fopen(fullfile(inpth,'output','output.gd'),'r') ;

    l = fgetl(fid) ;
    while ~contains(l,'#=PROGRAM')
        l = fgetl(fid) ;
    end
    lc = strsplit(l) ;
    tstf = find(~cellfun('isempty',(strfind(lc,'breseq')))) ;
    breseq_ver = lc{tstf+1} ;
    tstf = find(~cellfun('isempty',(strfind(lc,'revision')))) ;
    breseq_revision = lc{tstf+1} ;
    save (filename,"breseq_ver","breseq_revision")
end