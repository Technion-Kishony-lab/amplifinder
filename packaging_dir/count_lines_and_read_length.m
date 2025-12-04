function [num_lines,line_length,line_length_stats] = count_lines_and_read_length(fn)

line_length_stats = nan; 
if strcmp(fn(end-2:end),'.gz')
    l = evalc(['!zcat ' fn ' | tee >(wc -l) >(head -2 | tail -1 | wc -m) >/dev/null']);
    tstf = strfind(l,newline);
    n(1) = str2double(l(1:tstf(1)-1));
    n(2) = str2double(l(tstf(1)+1:tstf(2)-1));
    num_lines = max(n);
    line_length = min(n)-1;
else
    l = evalc(['!wc -l ' fn]);
    tstf = strfind(l,' ');
    num_lines = str2double(l(1:tstf(1)-1));
    l = evalc(['!head -2 ' fn ' | tail -1 | wc -m']);
    line_length = str2double(l)-1;
    v = [2:4:400];
    for i = 1:length(v)
        l = evalc(['!head -' num2str(v(i)) ' ' fn ' | tail -1 | wc -m']);
        line_length_stats(i) = str2double(l)-1; 
    end
end


