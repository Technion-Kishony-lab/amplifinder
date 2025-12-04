function [e,h] = parseBreseqCOV(iso,ref,inpth,outpth)

e = false ; 
h = nan;
ref = strsplit(ref,',') ;
missing = zeros(length(ref),1) ; 
for i = 1:length(ref)
   if ~exist(fullfile(outpth, [ref{i} '_COV.mat']),'file')
       missing(i) = 1 ; 
   end
end
if any(missing)
    tic
    %% Set up the Import Options and import the data
    opts = delimitedTextImportOptions("NumVariables", 8);
    
    % Specify range and delimiter
    opts.DataLines = [2, Inf];
    opts.Delimiter = "\t";
    
    % Specify column names and types
    opts.VariableNames = ["unique_top_cov", "unique_bot_cov", "Var3", "Var4", "Var5", "Var6", "Var7", "position"];
    opts.SelectedVariableNames = ["unique_top_cov", "unique_bot_cov", "position"];
    opts.VariableTypes = ["int32", "int32", "string", "string", "string", "string", "string", "double"];
    
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    
    % Specify variable properties
    opts = setvaropts(opts, ["Var3", "Var4", "Var5", "Var6", "Var7"], "WhitespaceRule", "preserve");
    opts = setvaropts(opts, ["Var3", "Var4", "Var5", "Var6", "Var7"], "EmptyFieldRule", "auto");
     
    for i = 1:length(ref)
           t = readtable(fullfile(inpth,'08_mutation_identification',[ref{i} '.coverage.tab']),opts) ; 
            h = height(t);
           if ~istable(t) || h<100
            e = true ; 
            %warning(['error reading coverage of isolate  ' iso{1}])
        else
            COV = t.unique_top_cov + t.unique_bot_cov ;
            T = toc ; 
            save (fullfile(outpth, [ref{i} '_COV.mat']),'COV','T')
        end
    end
end