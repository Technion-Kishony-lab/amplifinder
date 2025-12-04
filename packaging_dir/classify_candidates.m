function ISJC2 = classify_candidates(ISJC2)
global prms
%overall view: (same as junction2fasta)
%    ~~~>>>======>>>======>>>~~~
%(1)  ~~==
%(2)  ~~>>
%(3)           ==>>
%(4)            ====
%(5)              >>==
%(6)                      >>~~
%(7)                      ==~~

if height(ISJC2)<1
    return
end

JCT_PATTERN_NAMES = {
    % pattern multiple_copies name
    0 1 1 0 1 1 0  nan "flanked"                    % 1
    0 1 1 0 1 0 1  1 "hemi-flanked left"            % 2
    1 0 1 0 1 1 0  1 "hemi-flanked right"           % 3
    1 0 1 0 1 0 1  1 "unflanked"                    % 4
    0 1 0 0 1 0 1  0 "hemi-flanked left single"     % 5
    1 0 1 0 0 1 0  0 "hemi-flanked right single"    % 6
    1 0 0 0 0 0 1  0 "no IS"                        % 7
    0 1 0 0 0 1 0  0 "deletion"                     % 8
    };

JCT_PATTERN = cell2mat(JCT_PATTERN_NAMES(:, 1:7));
JCT_ISAMP = cell2mat(JCT_PATTERN_NAMES(:, 8));
JCT_NAMES = JCT_PATTERN_NAMES(:, 9);

ISO_ANC_PATTERN = {
    % iso anc left right append
    1 5 0 1 %"de novo right"
    1 6 1 0 %"de novo left"
    1 7 1 1 %"de novo left and right"
    2 5 0 0 %""
    2 7 1 0 %"de novo left"
    3 6 0 0 %""
    3 7 0 1 %de novo right
    4 7 0 0 %""
    };
ISO_ANC_PAIR = cell2mat(ISO_ANC_PATTERN(:, 1:2));
ISO_ANC_DENOVO = cell2mat(ISO_ANC_PATTERN(:, 3:4));

for i = 1:height(ISJC2)
    [name(i,1),iso_name(i,1),anc_name(i,1)] = classify_case(ISJC2(i,:));
end

ISJC2 = addvars(ISJC2,name,'NewVariableNames','event');
ISJC2 = addvars(ISJC2,iso_name,'NewVariableNames','isolate_architecture');
ISJC2 = addvars(ISJC2,anc_name,'NewVariableNames','ancestor_architecture');


    function [name,iso_name,anc_name] = classify_case(isjc2)
        [isok_iso, iso_vec] = get_jct_vec(isjc2.iso_jc_cov_left{1}, ...
            isjc2.iso_jc_cov_right{1}, isjc2.iso_jc_cov_green{1});
        [isok_anc, anc_vec] = get_jct_vec(isjc2.anc_jc_cov_left{1}, ...
            isjc2.anc_jc_cov_right{1}, isjc2.anc_jc_cov_green{1});

        if ~isok_iso || ~isok_anc
            low_coverage_near_jc = true;
        else
            low_coverage_near_jc = false;
        end

        pattern_iso = get_jct_pattern(iso_vec);
        pattern_anc = get_jct_pattern(anc_vec);

        if isempty(pattern_iso)
            name = "Unresolved";
            iso_name = "Unresolved";
        else
            name = JCT_NAMES{pattern_iso};
            iso_name = JCT_NAMES{pattern_iso};
        end

        if isempty(pattern_anc)
            name = "Unresolved";
            anc_name = "Unresolved";
        else
            anc_name = JCT_NAMES{pattern_anc};
        end

        if low_coverage_near_jc
            name = name + " (low coverage near junction)";
        end

        if ~isempty(pattern_anc) && ~isempty(pattern_iso)

            if pattern_iso==pattern_anc
                name = name + " (ancestral)";
                return
            end

            comb = find(all(ISO_ANC_PAIR==[pattern_iso,pattern_anc],2));
            assert(numel(comb)<=1)
            if isempty(comb)
                name = "Unresolved isolate-ancestor pair";
                return
            end
            if ISO_ANC_DENOVO(comb,1)
                name = name + " de novo left";
            end
            if ISO_ANC_DENOVO(comb,2)
                name = name + " de novo right";
            end
        end
    end

    function [is_ok, jct_vec] = get_jct_vec(left, right, green)
        is_ok = all(left >= prms.MIN_JCT_COV) && all(right >= prms.MIN_JCT_COV);
        jct_vec = green>=prms.MIN_JCT_COV;
    end

    function id = get_jct_pattern(jct_vec)
        id = find(all(jct_vec == JCT_PATTERN, 2));
        assert(numel(id)<=1)
    end

end