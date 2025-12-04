function n = count_lines(fn)

if strcmp(fn(end-3:end),'.gz')
    l = evalc(['zcat ' fn ' | wc -l']);
else
    l = evalc(['!wc -l ' fn]);
end


tstf = strfind(l,' ') ;
n = str2double(l(1:tstf(1)-1)) ;