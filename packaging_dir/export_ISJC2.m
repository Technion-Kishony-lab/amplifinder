function export_ISJC2(ISJC2,iso)

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
end
ISJC2 = addvars(ISJC2,ref','Before',"pos_Chr",'NewVariableNames',"Reference");
%output_ISJC2 = addvars(output_ISJC2,ref','After',"iso",'NewVariableNames',"Reference");
ISJC2 = addvars(ISJC2,ancestor','After',"amplicon_coverage_mode",'NewVariableNames',"Ancestor");
ISJC2 = addvars(ISJC2,IS','After',"dir_Chr",'NewVariableNames',"IS_element");
ISJC2 = addvars(ISJC2,DR','After',"IS_element",'NewVariableNames',"IS_direction");
ISJC2 = sortrows(ISJC2,"amplicon_coverage_mode",{'descend'});

ISJC2 = ISJC2(:,["iso","Reference","pos_Chr","dir_Chr","amplicon_length","IS_element","IS_direction","amplicon_coverage","amplicon_coverage_mode","Ancestor"]);
ISJC2 = renamevars(ISJC2,["pos_Chr","dir_Chr","amplicon_coverage","amplicon_coverage_mode"],...
                                        ["Positions_in_chromosome","Direction_in_chromosome","median_copy_number","mode_copy_number"]);

writetable(ISJC2,[iso.iso_outpath{1}, filesep, 'ISJC2.xlsx'])
z = (ISJC2.mode_copy_number>prms.COPY_NMBR_THRS | ISJC2.mode_copy_number<prms.DEL_COPY_NMBR_THRS) & ISJC2.amplicon_length>prms.FILTER_AMPLICON_LENGTH;
writetable(ISJC2(z,:),[iso.iso_outpath{1}, filesep, 'candidate_amplifications.xlsx'])
 
end