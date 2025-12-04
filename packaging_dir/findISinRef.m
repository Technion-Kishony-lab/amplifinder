function findISinRef(ref,ref_path)
% Find IS elements in genebank.
% Creates xxx_IS_loc.mat

IS_loc_file = [char(ref_path) filesep 'genbank' filesep ref '_IS_loc.mat'] ;
if exist(IS_loc_file,'file')
    return
end

out_raw_file = [char(ref_path) filesep 'genbank' filesep ref '_IS_loc_raw.txt'] ;
eval(['!grep -B1 "insertion sequence" ' char(ref_path) '/genbank/' ref '.gb > ' out_raw_file])
fraw = fopen(out_raw_file,'r') ;
IS_loc = table('Size', [0 7],'VariableTypes',{'double','string','string','double','double','logical','logical'}, ...
    'VariableNames',{'ID','IS_Name','IS_scaf','LocLeft','LocRight','Complement','Join'}) ;
k = 0;
while ~feof(fraw)
    l = fgetl(fraw) ;
    if k==0 && length(l)==1 && l==-1
        break
    end
    if contains(l,'--') || contains(l,'/') || contains(l,'=')
        continue
    end
    k = k + 1;
    scaf = ref ;
    fstmp1 = strfind(l,'.') ;
    re = regexp(l,'[0-9]') ;
    posL = str2num(l(re(1):fstmp1(1)-1)) ;
    posR = str2num(l(fstmp1(end)+1:re(end))) ;
    comp = contains(l,'complement');
    join = contains(l,'join');
    l = fgetl(fraw) ;
    fstmp1 = strfind(l,':');
    if ~isempty(fstmp1)
        fstmp2 = strfind(l(fstmp1+1:end),'"');
        name = l(fstmp1+1:fstmp1+fstmp2-1) ;
    else
        name = 'unknown' ;
    end
    IS_loc = [IS_loc; {k,name,scaf,posL,posR,comp,join}] ;
    fgetl(fraw) ;
end
save(IS_loc_file, 'IS_loc')

end
