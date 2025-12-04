function ISJC2 = curate_candidate_amplicons(iso,anc_outpath)

% ISJC: IS-associated junction. Junctions can be based on expected IS and
% chromosome junctions according to reference genome, or based on breseq
% junction analysis.
% ISJC2: table of all pairs of ISJC potentially of the same IS element and facing in opposite directions.
% JC fields: parsed from Breseq
% ISJC fields: JC fields + IS_side, IS_diff, switched
% ISJC2 fields: iso, JC_num, pos_Chr, pos_IS, dir_Chr, dir_IS, ref_IS,
% read_count_offset, IS_distance, breseq_reject, scaf_Chr, IS,
% span_origin, amplicon_length, complementary_length, amplicon_coverage,
% amplicon_coverage_mode, iso_outpath, anc_outpath

global prms

%% Only curating ISJC2 if .mat file does not already exist
iso_outpath = iso.iso_outpath{1};
filename = [iso_outpath filesep 'ISJC2.mat'];
if exist(filename,'file')
    load (filename,"ISJC2")
    return
end

%% Constructing ISJC based on previously generated JC and JC_refIS
load(fullfile(iso_outpath, 'MUT.mat'), 'JC') ;
load(fullfile(iso_outpath, 'ref_props.mat'), 'ref_props');

% JC
JC = addvars(JC,zeros(height(JC),1),'After','dir2','NewVariableNames','refIS');
% Decision to remove specific JCs based on breseq criterion of rejection.
remove_JC_breseq_reject = prms.REMOVE_JC_BRESEQ_REJECT;
if remove_JC_breseq_reject
    JC = JC(ismissing(JC.reject), :); % remove BRESEQ-rejected JCs
end
JC = replace_scaf_name_with_index(JC, ref_props.name);

% JC_refIS
load(fullfile(iso_outpath, 'JC_refIS'), 'JC_refIS');

% JC + JC_refIS
allJC = [JC; JC_refIS]; % concatenate de novo JCs and JCs generated from reference
% assigning ISs
ISJC = assign_potential_ISs(allJC,iso.ref_path(1),ref_props.name,iso.isfinder,prms.MAX_DIST_TO_IS);

%% curating ISJC2 from ISJC, also using ref_props and COV
ISJC2 = combine_ISJC_pairs(ISJC);
if height(ISJC2)>0
    ISJC2 = calculate_amplicon_length(ISJC2, ref_props);
end

if height(ISJC2)>0
    [amplicon_copy_number_dist,ISJC2] = calc_coverage_ISJC2(ISJC2, ref_props, iso_outpath, anc_outpath);
else
    amplicon_copy_number_dist = 0;
end
X1 = prms.NCP_LIMIT1;
X2 = prms.NCP_LIMIT2;

% add columns: iso, iso_outpath, anc_outpath, ref_path
h = height(ISJC2);
ISJC2 = [ISJC2, iso(ones(h,1), {'iso','iso_outpath','ref_path'})];
ISJC2 = movevars(ISJC2,"iso",'Before',"JC_num"); %iso column becomes first one.
ISJC2 = addvars(ISJC2, repmat(anc_outpath,h,1), 'NewVariableNames','anc_outpath');

ISJC2 = [ISJC2, rowfun(@get_juction_name, ISJC2(:,{'scaf_Chr','IS', 'dir_Chr', 'pos_Chr'}), ...
    'OutputVariableNames','analysis_directory')];

save([iso_outpath filesep 'ISJC.mat'],'ISJC','anc_outpath')
save(filename,'ISJC2','anc_outpath','amplicon_copy_number_dist','X1','X2');

logger(fprintf('Generated ISJC2 for %s\n', iso.iso{1}),'MESSAGE')

    function name = get_juction_name(scaf_Chr, IS, dir_Chr, pos_Chr)
    chr_name = ref_props.name(scaf_Chr);
    IS = IS{1}(1);
    IS_name = num2str(IS.id);
    IS_orientation = orientation_to_char(IS.orientation);
    dir_Chr = orientation_to_char(dir_Chr);
    chr_jcn_name = sprintf('chr_%s_%d%s_%d%s',chr_name, pos_Chr(1), dir_Chr(1), pos_Chr(2), dir_Chr(2));
    IS_jcn_name = sprintf('IS_%s%s', IS_name, IS_orientation);
    name = string([chr_jcn_name '_' IS_jcn_name]);
    end

end

function c = orientation_to_char(orientation)
dir_str = 'RUF';  % Forward, Unknown, Reverese
c = dir_str(orientation + 2);
end