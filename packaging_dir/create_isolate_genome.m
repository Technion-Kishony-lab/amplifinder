function create_isolate_genome(ref_path,ref_names,output_path)

filename = strcat(output_path, filesep,"ref_props.mat");
if ~exist(filename,'file')
    files = cellfun(@(ref_name) [char(ref_path) '/genbank' filesep ref_name 'ref.mat'], ref_names,'UniformOutput',false);
    ref_props = concatenateTableFiles(files);
    save (filename,"ref_props")
end

end
