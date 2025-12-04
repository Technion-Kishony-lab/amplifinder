function prms = readConfig()
    
    global workingDir baseDir
    
    config_file_bame = 'config.txt';
    % Define the filename
    filename = fullfile(workingDir, config_file_bame);
    if exist(filename, "file")
        fprintf('Overrding default config with local config.txt')
    else
        filename = fullfile(baseDir, config_file_bame);
    end
    
    % Check if the file exists
    if ~isfile(filename)
        fprintf('File "config.txt" not found.');
    end
    
    % read parameters (resulting in lowercase
    p = ini2struct(filename);
    flds = fields(p);
    prms = struct(); 
    for i = 1:length(flds)
        fieldname = flds{i};
        if ~endsWith(fieldname, 'path')
            prms.(upper(fieldname)) = str2num(p.(fieldname));
        else
            prms.(upper(fieldname)) = p.(fieldname); 
        end
    end

    % Display the parsed configuration
    fprintf('Configuration loaded from "%s":\n', filename);
    disp(prms);
end
