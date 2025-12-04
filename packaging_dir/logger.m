function logger(msg, type, should_print, should_log)
if nargin<2
    type = 'MESSAGE';
end
if nargin<3
    should_print = true;
end
if nargin<4
    should_log = true;
end

isError = strcmp(type, 'ERROR');
global prms logDir
filepath = fullfile(logDir, prms.LOG_PATH);

if isstring(msg) || ischar(msg)
    full_msg = sprintf('%s: %s | %s\n', type,msg,char(datetime('now')));
else
    full_msg = evalc('disp(msg)');
end

if should_print && ~isError
    fprintf(full_msg)
end

if should_log
    fid = fopen(filepath, "a");
    fprintf(fid, full_msg);
end
fclose(fid);

if isError
    error(full_msg)
end

end