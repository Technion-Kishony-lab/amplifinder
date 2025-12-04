function JC = replace_scaf_name_with_index(JC, ref_names)

[~,JC.jscaf1] = ismember(JC.scaf1, ref_names);
[~,JC.jscaf2] = ismember(JC.scaf2, ref_names);

JC.scaf1 = [];
JC.scaf2 = [];

end
