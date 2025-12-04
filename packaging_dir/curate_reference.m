function iso = curate_reference(iso,NCBI)
% assumes all references are in the same directory (all iso.ref_path are identical)

global prms
iso = addvars(iso,strings(height(iso),1),'NewVariableNames','ref','After','input_ref');
iso = addvars(iso,zeros(height(iso),1),'NewVariableNames','number_IS_genbank','After','ref');
iso = addvars(iso,zeros(height(iso),1),'NewVariableNames','number_IS_ISfinder','After','number_IS_genbank');
input_ref_names = cellfun(@(x)split(x,','), iso.input_ref, 'UniformOutput',false);
unique_input_ref_names = unique(cat(1, input_ref_names{:}));
if NCBI
    ref_names = arrayfun(@(x, y) get_reference(x{1}, y{1}), unique_input_ref_names, repmat(cellstr(iso.ref_path(1)),size(unique_input_ref_names,1),size(unique_input_ref_names,2)), 'UniformOutput', false);
    unique_ref_names = unique(ref_names);
else
    ref_names = unique_input_ref_names;
    unique_ref_names = ref_names;
end

cellfun(@(ref_name)findISinRef(ref_name,iso.ref_path(1)),ref_names)
cellfun(@(ref_name)ISfinder(ref_name,iso.ref_path(1)),ref_names)

cellfun(@(ref)create_IS_end_seqs(ref,iso.ref_path(1),prms.LENGTH_SEQ_INTO_IS,prms.MAX_DIST_TO_IS), unique_ref_names);
fprintf('\n')
for i = 1:height(iso)
    [~,j] = ismember(input_ref_names{i}, unique_input_ref_names);
    iso.ref(i) = strjoin(ref_names(j),',');
    create_isolate_genome(iso.ref_path(1),unique_input_ref_names(j), iso.iso_outpath(i))
    [iso.number_IS_genbank(i),iso.number_IS_ISfinder(i)] = create_isolate_ISloc(iso.ref_path(1),ref_names(j), iso.iso_outpath(i), iso.isfinder(i));
    create_JC_of_reference_IS(iso.iso_outpath(i));
end

if (iso.number_IS_genbank(1)==0 && ~iso.isfinder(1)) || (iso.number_IS_ISfinder(1)==0 && iso.isfinder(1))
    logger('no IS elements identified in genome. Exiting...','ERROR')
end