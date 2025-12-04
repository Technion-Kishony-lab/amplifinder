function IS_loc = blast2ISloc(blast_results,headers,len,critical_coverage)
if nargin<4
    critical_coverage = 0.9;
end

t = blast_results;

if height(t) > 0
    if height(t) > 1
        [u,it,iu] = unique([t.qstart,t.qend],"rows",'first') ;
        t = t(it,:) ;
    end
    [~,fm] = ismember(t.subject,headers) ;
    len = len(fm) ;
    t = addvars(t,len,'NewVariableNames','subject_length','After','length') ;
    t = sortrows(t,{'query','qstart','bitscore'},'ascend') ;
    qstart = t.qstart(2:end) ;
    qend = t.qend(1:end-1) ;
    d = [1 ; qstart - qend ] ;
    f_pos = [find(d>0);length(d)+1] ;
    k = false(height(t),1) ;
    for i = 1:length(f_pos)-1
        o = f_pos(i):f_pos(i+1)-1 ;
        [m,z] = min(abs(t.qend(o)-t.qstart(o)-t.subject_length(o))) ;
        if abs(t.sstart(o(z)) - t.send(o(z)))/t.subject_length(o(z)) > critical_coverage
            k(o(z)) = true ;
        end
    end
    t = t(k,:) ;
end

%generating IS_loc
IS_loc = t(:, {'subject', 'query', 'qstart', 'qend'});
IS_loc.Properties.VariableNames = {'IS_Name', 'IS_scaf', 'LocLeft', 'LocRight'};
IS_loc.Complement = false(height(t),1) ;
IS_loc.Join = false(height(t),1) ;
IS_loc = addvars(IS_loc, (1:height(IS_loc))', 'NewVariableNames','ID', 'Before','IS_Name');
