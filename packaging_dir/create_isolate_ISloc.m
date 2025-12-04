function [n_genbank,n_ISfinder] = create_isolate_ISloc(ref_path,ref_names,output_path,isfinder)

filepath = strcat(output_path, filesep, "num_IS_loc.mat");

source = {'genbank','ISfinder'};
    if isfinder == 0
        nsource = 1;
    else
        nsource = 2;
    end

if exist(strcat(output_path, filesep, "IS_loc.mat"),'file')
    load (strcat(output_path, filesep, "IS_loc.mat"),'IS_loc_source')
    if ~exist('IS_loc_source','var') || ~strcmp(IS_loc_source,source(nsource))
        delete(strcat(output_path, filesep, "IS_loc.mat"))
        delete(strcat(output_path, filesep, "num_IS_loc.mat"))
        fprintf('recalculating IS_loc in %s\n',output_path);
    end
end

if ~exist(filepath,'file')
    for i = 1:length(source)
        files = cellfun(@(ref_name) [char(ref_path) filesep source{i} filesep ref_name '_IS_loc.mat'], ref_names,'UniformOutput',false);
        IS_loc = concatenateTableFiles(files);
        if height(IS_loc)==0
            fprintf('IS_loc %s empty for isolate in path %s\n',source{i}, output_path)
        else
            % replace scaffold name with index:
            [~, IS_loc.IS_jscaf] = ismember(IS_loc.IS_scaf, ref_names);
            IS_loc.IS_scaf = [];
        end
        save(strcat(output_path, filesep, source{i}, "_IS_loc.mat"), 'IS_loc');
        if i==1
            n_genbank = height(IS_loc);
        else
            n_ISfinder = height(IS_loc);
        end
        if nsource == i
            IS_loc_source = source(nsource);
            save(strcat(output_path, filesep, "IS_loc.mat"), "IS_loc","IS_loc_source");
        end
    end
    save(filepath, "n_genbank","n_ISfinder");
else
    load(filepath, "n_genbank","n_ISfinder");
end
end
