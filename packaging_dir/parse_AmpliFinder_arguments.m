function [iso,NCBI,steps] = parse_AmpliFinder_arguments(ARGS)

% Default values for arguments
iso = table();
iso.iso(1)="isolate";                           %default for -iso_name
iso.iso(2)="ancestor"; iso.anc(1)="ancestor";   %default for -anc_nam
iso.anc(2)="";                                  %second isoalte is ancestor and has no defined ancestor
iso.ref_path(1)="genomesDB"; iso.ref_path(2)="genomesDB"; %default for -ref_path
iso.isfinder(1)=false; iso.isfinder(2)=false;   %default for -ISfinder 
NCBI=true;                                      %default for -NCBI
iso.breseq_path(1) = ""; iso.breseq_path(2) = ""; 
steps = "all";

% List of required flags
requiredFlags = {'-iso_path','-ref_name'};
providedFlags = {};  % Keep track of provided flags

% Argument parsing
i = 1;
while i <= length(ARGS)
    switch ARGS{i}
        case '-iso_path'
            % Check if a value follows the flag
            if i + 1 <= length(ARGS)
                iso.fastq_path(1) = string(ARGS{i + 1});
                providedFlags{end+1} = '-iso_path';  % Record flag as provided
                i = i + 2; % Skip the value after the flag
            else
                logger('Missing value for argument "-iso_path".','ERROR');
            end
        case '-iso_breseq_path'
            % Check if a value follows the flag
            if i + 1 <= length(ARGS)
                iso.breseq_path(1) = string(ARGS{i + 1});
                providedFlags{end+1} = '-iso_breseq_path';  % Record flag as provided
                i = i + 2; % Skip the value after the flag
            else
                logger('Missing value for argument "-iso_breseq_path".','ERROR');
            end
        case '-iso_name'
            % Check if a value follows the flag
            if i + 1 <= length(ARGS)
                iso.iso(1) = string(ARGS{i + 1});
                providedFlags{end+1} = '-iso_name';  % Record flag as provided
                i = i + 2; % Skip the value after the flag
            else
                logger('Missing value for argument "-iso_name".','ERROR');
            end
        case '-anc_path'
            % Check if a value follows the flag
            if i + 1 <= length(ARGS)
                iso.fastq_path(2) = string(ARGS{i + 1});
                providedFlags{end+1} = '-anc_path';  % Record flag as provided
                i = i + 2; % Skip the value after the flag
            else
                logger('Missing value for argument "-anc_path".','ERROR');
            end
        case '-anc_breseq_path'
            % Check if a value follows the flag
            if i + 1 <= length(ARGS)
                iso.breseq_path(2) = string(ARGS{i + 1});
                providedFlags{end+1} = '-anc_breseq_path';  % Record flag as provided
                i = i + 2; % Skip the value after the flag
            else
                logger('Missing value for argument "-anc_breseq_path".','ERROR');
            end
        case '-anc_name'
            % Check if a value follows the flag
            if i + 1 <= length(ARGS)
                iso.iso(2) = string(ARGS{i + 1});
                iso.anc(1) = string(ARGS{i + 1});
                providedFlags{end+1} = '-anc_name';  % Record flag as provided
                i = i + 2; % Skip the value after the flag
            else
                logger('Missing value for argument "-anc_name".','ERROR');
            end
        case '-ref_path'
            % Check if a value follows the flag
            if i + 1 <= length(ARGS)
                iso.ref_path(1) = string(ARGS{i + 1});
                iso.ref_path(2) = string(ARGS{i + 1});
                providedFlags{end+1} = '-ref_path';  % Record flag as provided
                i = i + 2; % Skip the value after the flag
            else
                logger('Missing value for argument "-ref_path".','ERROR');
            end
        case '-ref_name'
            % Check if a value follows the flag
            if i + 1 <= length(ARGS)
                iso.input_ref(1) = string(ARGS{i + 1});
                iso.input_ref(2) = string(ARGS{i + 1});
                providedFlags{end+1} = '-ref_name';  % Record flag as provided
                i = i + 2; % Skip the value after the flag
            else
                logger('Missing value for argument "-ref_name".','ERROR');
            end
        case '-steps'
            % Check if a value follows the flag
            if i + 1 <= length(ARGS)
                steps = string(ARGS{i + 1});
                providedFlags{end+1} = '-steps';  % Record flag as provided
                i = i + 2; % Skip the value after the flag
            else
                logger('Missing value for argument "-steps".','ERROR');
            end
        case '-NCBI'
            % Check if a value follows the flag
            if i + 1 <= length(ARGS)
                if strcmp(ARGS{i + 1},"true") | strcmp(ARGS{i + 1},"TRUE")
                    NCBI = logical(true);
                elseif strcmp(ARGS{i + 1},"false") | strcmp(ARGS{i + 1},"FALSE")
                    NCBI = logical(false);
                else
                    logger('Argument "-NCBI" expects "true" or "false" value.','ERROR');
                end
                providedFlags{end+1} = '-NCBI';  % Record flag as provided
                i = i + 2; % Skip the value after the flag
            else
                logger('Missing value for argument "-NCBI".','ERROR');
            end
        case '-ISfinder'
            % Check if a value follows the flag
            if i + 1 <= length(ARGS)
                if strcmp(ARGS{i + 1},"true") | strcmp(ARGS{i + 1},"TRUE")
                    iso.isfinder(1) = logical(true);
                    iso.isfinder(2) = logical(true);
                elseif strcmp(ARGS{i + 1},"false") | strcmp(ARGS{i + 1},"FALSE")
                    iso.isfinder(1) = logical(false);
                    iso.isfinder(2) = logical(false);
                else
                    error('Argument "-ISfinder" expects "true" or "false" value.');
                end
                providedFlags{end+1} = '-ISfinder';  % Record flag as provided
                i = i + 2; % Skip the value after the flag
            else
                logger('Missing value for argument "-ISfinder".','ERROR');
            end
        otherwise
            logger(fprintf('Unknown argument "%s".', ARGS{i}),'ERROR');
    end
end

% Defining outpath and breseqpath
for i = 1:height(iso)
    iso.iso_outpath(i) = strcat("output",filesep,iso.iso(i));
    if iso.breseq_path(i) == ""
        iso.breseq_path(i) = strcat('BRESEQ',filesep, iso.iso(i));
    end
end

% Validate that all required flags are provided
missingFlags = setdiff(requiredFlags, providedFlags);
if ~isempty(missingFlags)
    logger(fprintf('Missing required arguments: %s', strjoin(missingFlags, ', ')),'ERROR');
end

% if no -anc_path, remove ancestor line from iso table
if ~ismember('-anc_path',providedFlags) && ~ismember('-anc_name',providedFlags)
    iso = iso(1,:);
    logger('No ancestor assigned; using non-normalized coverage analysis','WARNING');
end

% genomesDB path is reserved for downloading from NCBI
if strcmp(iso.ref_path(1),"genomesDB") && ~NCBI
    logger('genomesDB directory is reserved for ncbi references. Choose a different directory.','ERROR')
end

% ISfinder is required for non-NCBI (local) genomes
if ~iso.isfinder(1) && ~NCBI
    logger('ISfinder is required for local (non-ncbi) genomes. Use ''-ISfinder true''.','ERROR' )
end
end