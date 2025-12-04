function ISJC2 = classify_ISJC2(ISJC2)
% classify ISJC2 as (1) internal-only (no flanking IS elements), (2) hemi-flanked, (3) flanked, based on 
% the contributaion of their ISJCs to single-locus IS elements (either transposition or reference IS).  

% identify canonical and non-caninical transpositions
%
% legend:
%      >>> IS
%   ====== amplified region (cassette)
%     ~~~~ flanking chromosome
%       -- IS-associated junction (ISJC)
%
%
% transposition:
% ~~~AB~~~                       reference
% ~~~A>>>B~~~                    isolate
%    -- --                       focal junction-pair
%
% ancestral IS: 
% ~~~~>>>~~~~                    reference
% ~~~~>>>~~~~                    isolates
%    -- --                       focal junction-pair
% 
% unflanked structure:    both ISJCs are not contributing to a single-locus IS                            
% ~~~Aa======bB~~~               reference
% ~~~Aa======b>>>a======bB~~~    isolate
%            -- --               focal junction-pair
%
% hemi-flanked structure:
% ~~~Aa======bB~~~               reference
% ~~~A>>>a======b>>>B=======~~~~ isolate
%       --      --               focal junction-pair
%    -- --                       one ISJC also contributes to a transposition
%
% hemi-flanked with ancestral IS:
% ~~~Aa======b>>>B~~~            reference
% ~~~Aa======b>>>a======b>>>B~~~ isolate
%               --      --       focal junction-pair
%                       -- --    one ISJC also contributes to a reference IS
%                   
% flanked structure with transpositions/reference IS:
% ~~~Aa======b>>>B~~~            reference
% ~~~A>>>a======b>>>B~~~         isolate
%       --      --               focal junction-pair
%    -- --                       transposition junction-pair 1
%               -- --            reference IS junction-pair 2
% 
% position of cassette relative to origin
% legend:
% L - position of left chromosomal position of junction-pair
% R - position of right chromosomal position of junction-pair
%
% not spanning origin:
% ~~~~L======R~~~~o~~~~         reference
%
% spanning origin:
% ~~~~R====o===L~~~~            reference
%

% maximal distance between insertion sites of two sides of IS to be 
% considerted as transposition:

global prms

TRANSPOSITION_THRESHOLD = prms.MIN_AMPLICON_LENGTH;  %30

h = height(ISJC2);

is_reference = all(ISJC2.refIS>0, 2) & diff(ISJC2.refIS,[],2)==0; %a transposon in its reference locus.
is_hemi_reference = sum(ISJC2.refIS>0, 2)==1; %one of the pair of junctions belongs to a transposon in its reference locus.
is_denovo = ~is_reference & ~is_hemi_reference;
is_transposition =  is_denovo & ISJC2.amplicon_length<TRANSPOSITION_THRESHOLD; % transposition made of two closely located junctions. 

is_single_locus = is_transposition | is_reference; %all transposons, denovo or reference. 

% for each IS-Chr junction, we check if it ever appears in a "normal" IS 
% (namly is it ever a part of a junction pair that is either reference or transposition). 
x = directed_IS;
empty_directed_IS = x(1,[]);
ISs_shared_with_single_locus_pair = repmat({empty_directed_IS}, h, 2);
ISs_shared_with_both_single_locus_pairs = repmat({empty_directed_IS}, h,1);
pos_of_paired_single_locus_junction = zeros(h,2,'int32');
dis_of_paired_single_locus_junction = cell(h,2);
IS_of_paired_single_locus_junction = cell(h,2);
multiple_single_locus = false(h,1);
for i = 1:h
    for j = 1:2
        JC_num_ij = ISJC2.JC_num(i,j); %reminder: 0 for reference, integer for breseq id
        ref_ij = ISJC2.refIS(i,j); %reminder: 0 for non-reference, integer for the index of the reference transposon.
        dir_IS_ij = ISJC2.dir_IS(i,j);
        %is_same_as_focal_jcn: array, size table height X 2, indicating where the 'focal' junction appears 
        if JC_num_ij>0 % non-reference JC
            is_same_as_focal_jcn = ISJC2.JC_num==JC_num_ij; 
        else
            is_same_as_focal_jcn = ISJC2.refIS==ref_ij & ISJC2.dir_IS==dir_IS_ij;
        end
        is_same_as_focal_jcn_and_single_locus_pair = is_single_locus & is_same_as_focal_jcn; %returns an array
        is_same_as_focal_jcn_and_single_locus_pair(i,:) = false;  % remove self
        [ii,jj] = find(is_same_as_focal_jcn_and_single_locus_pair);
        %assert(numel(ii)<=1)  % a junction can only appear in one single locus pair
        if numel(ii)>1
            multiple_single_locus(i) = true; % multiple potential single locus for the same junctions are treated as a "bug" or an artifact. Ignored. 
        elseif numel(ii)==1
            jj_rev = 3-jj; % 1->2, 2->1

            IS_of_paired_single_locus_junction(i,j) = ISJC2.IS(ii); % a list of the IS elements of the single locus
            pos_of_paired_single_locus_junction(i,j) = ISJC2.pos_Chr(ii,jj_rev); %list
            dis_of_paired_single_locus_junction(i,j) = ISJC2.IS_distance(ii,jj_rev); %list

            % find matching IS (same ID and same direction)
            [matched_ISs,ia,ib] = intersect(ISJC2.IS{i}, ISJC2.IS{ii});
            ISs_shared_with_single_locus_pair{i,j} = matched_ISs;
        end
    end
    [ISs_shared_with_both_single_locus_pairs{i},ia,ib] = intersect(ISs_shared_with_single_locus_pair{i,:}); % IS elements shared by focal ISJC2 and the bordering IS elements if exist.
end

is_ISIDs_shared_with_both_single_locus_pairs = ~cellfun(@isempty, ISs_shared_with_both_single_locus_pairs);
is_ISIDs_shared_with_single_locus_pair = ~cellfun(@isempty, ISs_shared_with_single_locus_pair);
%is_ISIDs_shared_with_single_locus_pair = ~cellfun(@isempty, IS_of_paired_single_locus_junction);

pos_of_paired_single_locus_junction(~is_ISIDs_shared_with_single_locus_pair) = 0;

ISJC2.pos_of_paired_single_locus_junction = pos_of_paired_single_locus_junction;
ISJC2.dis_of_paired_single_locus_junction = dis_of_paired_single_locus_junction;
ISJC2.IS_of_paired_single_locus_junction = IS_of_paired_single_locus_junction;

event = repmat("unresolved", h,1);
event(multiple_single_locus) = "multiple single locus";
event(is_reference & ~multiple_single_locus) = "reference"; 
event(is_transposition & ~multiple_single_locus) = "transposition";

is_not_flanking = sum(is_ISIDs_shared_with_single_locus_pair,2)==0 & ~is_single_locus;
event(is_not_flanking & ~multiple_single_locus) = "unflanked";

is_hemi_flanking = sum(is_ISIDs_shared_with_single_locus_pair,2)==1 & ~is_single_locus;
paired_side = ~cellfun(@isempty,IS_of_paired_single_locus_junction);
is_hemi_flanking_left = is_hemi_flanking & ((paired_side(:,1) & ~ISJC2.span_origin) | (paired_side(:,2) & ISJC2.span_origin)); 
is_hemi_flanking_right = is_hemi_flanking & ((paired_side(:,2) & ~ISJC2.span_origin) | (paired_side(:,1) & ISJC2.span_origin));
event(is_hemi_flanking_left & ~multiple_single_locus) = "hemi-flanked left";
event(is_hemi_flanking_right & ~multiple_single_locus) = "hemi-flanked right";

is_flanking = is_ISIDs_shared_with_both_single_locus_pairs & ~is_single_locus;
event(is_flanking & ~multiple_single_locus) = "flanked";

% shared IS:
shared_directed_IS = cell(h,1);

% flanked
shared_directed_IS(is_flanking) = ISs_shared_with_both_single_locus_pairs(is_flanking);

% internal-only
shared_directed_IS(is_not_flanking) = ISJC2.IS(is_not_flanking);

% hemi-flanked
ok = is_ISIDs_shared_with_single_locus_pair & is_hemi_flanking; 

ff = find(any(ok,2));
for i = 1:length(ff)
    if ok(ff(i),1)
        shared_directed_IS(ff(i)) = ISs_shared_with_single_locus_pair(ff(i),1);
    else
        shared_directed_IS(ff(i)) = ISs_shared_with_single_locus_pair(ff(i),2);
    end
end

ISJC2.raw_event = event;
ISJC2.shared_IS = shared_directed_IS;

% chosen IS:
is_shared_IS = ~cellfun(@isempty, shared_directed_IS);
chosen_IS = repmat({directed_IS}, h, 1);
chosen_IS(is_shared_IS) = cellfun(@(x){x(1)}, shared_directed_IS(is_shared_IS));
chosen_IS(~is_shared_IS) = cellfun(@(x){x(1)}, ISJC2.IS(~is_shared_IS)); 
ISJC2.chosen_IS = chosen_IS;
end

