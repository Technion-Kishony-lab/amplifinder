function export_filteredISJC2(ISJC2,iso)

fn = [iso.iso_outpath{1}, filesep, 'classified_amplifications.xlsx'];
if height(ISJC2)<1
    save ([iso.iso_outpath{1}, filesep, 'classified_amplifications.mat'],'ISJC2')
    writetable(ISJC2,fn)
    return
end

global prms
load (char(strcat(iso.iso_outpath(1),filesep,"IS_loc")),"IS_loc")

clear ref ancestor
for i = 1:height(ISJC2)
    temp_ref = strsplit(iso.ref(1),',');
    ref(i) = temp_ref(ISJC2.scaf_Chr(i)); 
    ancestor(i) = iso.iso(2); 
    R = extractBefore([ISJC2.IS{i}.name],"R"); % IDs of reverse facing IS elements
    F = extractBefore([ISJC2.IS{i}.name],"F"); % IDs of forward facing IS elements
    NN = [F(~ismissing(F)),R(~ismissing(R))]; 
    DR(i) = strjoin([repmat("F",1,nnz(~ismissing(F))),repmat("R",1,nnz(~ismissing(R)))],','); 
    IS(i) = strjoin(IS_loc.IS_Name(str2double(NN),:),',');
    chosen_IS_explicit(i) = IS_loc.IS_Name(ISJC2.chosen_IS{i}.id);
end
ISJC2 = addvars(ISJC2,ref','Before',"pos_Chr",'NewVariableNames',"Reference");
%output_ISJC2 = addvars(output_ISJC2,ref','After',"iso",'NewVariableNames',"Reference");
ISJC2 = addvars(ISJC2,ancestor','After',"amplicon_coverage_mode",'NewVariableNames',"Ancestor");
ISJC2 = addvars(ISJC2,IS','After',"dir_Chr",'NewVariableNames',"IS_element");
ISJC2 = addvars(ISJC2,DR','After',"IS_element",'NewVariableNames',"IS_direction");
ISJC2 = addvars(ISJC2,chosen_IS_explicit','After',"IS_direction",'NewVariableNAmes',"chosen_IS_explicit");
ISJC2 = sortrows(ISJC2,"amplicon_coverage_mode",{'descend'});

ISJC2 = ISJC2(:,["iso","Reference","pos_Chr","dir_Chr","amplicon_length","IS_element","IS_direction","chosen_IS_explicit"...
    "amplicon_coverage","amplicon_coverage_mode","raw_event","Ancestor","analysis_directory","iso_jc_cov_green","iso_jc_cov_left","iso_jc_cov_right",...
    "anc_jc_cov_green","anc_jc_cov_left","anc_jc_cov_right","event","isolate_architecture","ancestor_architecture"]);
ISJC2 = renamevars(ISJC2,["pos_Chr","dir_Chr","amplicon_coverage","amplicon_coverage_mode","iso_jc_cov_green","anc_jc_cov_green","chosen_IS_explicit"],...
                         ["Positions_in_chromosome","Direction_in_chromosome","median_copy_number","mode_copy_number","iso_jc_cov_trans","anc_jc_cov_trans","sequence_IS"]);

save ([iso.iso_outpath{1}, filesep, 'classified_amplifications.mat'],'ISJC2')

%converting vectors to strings to avoid multiple columns in xlsx
converted_fields = ["Positions_in_chromosome","Direction_in_chromosome","iso_jc_cov_trans","iso_jc_cov_left","iso_jc_cov_right","anc_jc_cov_trans","anc_jc_cov_left","anc_jc_cov_right"];

for i = 1:length(converted_fields)
    if ~iscell(ISJC2.(converted_fields(i)))
        ISJC2.(converted_fields(i)) = mat2cell(ISJC2.(converted_fields(i)),ones(1,height(ISJC2)),2);
    end
    ISJC2.(converted_fields(i))= cellfun(@(x) strjoin(string(x), ','), ISJC2.(converted_fields(i)), 'UniformOutput', false);
end

writetable(ISJC2,fn)
 
end