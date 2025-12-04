function ISJC2 = combine_ISJC_pairs(ISJC)
% Create pairs of ISJC (junctions whose side "1" is IS) 
% if the following 3 conditions apply:
% (a) the two junctions are on the same chromosome facing opposing directions
% (b) the junctions match different sides of the same IS (or of multiple IS).


% (a1) chromosomal sides should be on same genome
same_scaf = ISJC.jscaf2 == ISJC.jscaf2';

% (a2) opposing chromosomal sides:
opposing_chromosome_dirs = ISJC.dir2 ~= ISJC.dir2';

all_ok = same_scaf & opposing_chromosome_dirs;

% Create empty ISJC2 table (table of pairs of JCs):
ISJC2 = table;

% Run over all pairs:
k = 0 ;
n = height(ISJC);
breseq_reject = ~ismissing(ISJC.reject) & ISJC.num~=0;
for i = 1:n-1
    for j = i+1:n
        %assert(ISJC.pos2(i) <= ISJC.pos2(j)) % i is used for L and j is used for R
        if ~all_ok(i,j)
            continue
        end
        
        ID_i = [ISJC.IS_side{i}.id];
        side_i = [ISJC.IS_side{i}.orientation]; % oreintation here means side of the IS (-1, 1)
        ID_j = [ISJC.IS_side{j}.id];
        side_j = [ISJC.IS_side{j}.orientation];

        % (b) find same IS different sides
        match = ID_i==ID_j' & xor(side_i==1,side_j'==1);
        
        % skip pair if no IS matches the two junctions 
        if ~any(match, [1,2])
            continue 
        end

        % span origin (assuming scaf circularity, which we check later):
        span_origin = ISJC.dir2(i)==-1;

        % find all pairs of matching IS sides:
        [jj,ii] = find(match);
        
        % skip pair if a reference junction is involved, and the list of matching ISs does not include the reference
        % however, IS_id is not limited to the reference IS. This limiting may result in false negative. e.g. a flanked amplification 
        % formed by a transposition of IS X and a reference Y.    
        ID_refs = unique([uint8(ISJC.refIS(i)),uint8(ISJC.refIS(j))]);
        ID_refs = ID_refs(ID_refs~=0);
        if ~isempty(ID_refs) && ~all(ismember(ID_refs,ID_i(ii)))
            continue
        end
        IS_id = ID_i(ii); % this equals ID_j(jj)

        % orientation of IS:     
        % 1: same as in reference. -1: opposite than reference. 0: match both directions. 
        %
        % if the IS is in its reference orientation,
        % then the left side of the amplicon should connect to the right
        % side of the IS:
        orientation = side_i(ii);

        % If spanning origin, then the orientation is inverted:
        orientation = orientation * ISJC.dir2(i); % spanning origin

        % Orientation
        % for ~span_origin:
        % ~~~~~1~L>>R~2~~~~  +1
        % ~~~~~1~R>>L~2~~~~  -1
        % ~~~~~1~B>>B~2~~~~   0
        %
        % for span__origin:
        % ~~~~~2~L>>R~1~~~~  -1
        % ~~~~~2~R>>L~1~~~~  +1
        % ~~~~~2~B>>B~1~~~~   0
        %

        % create a new line with joint columns:
        k = k + 1 ;
        ij = [i,j];
        m = [1 2];

        % Store unique properties of the two junctions:
        ISJC2.JC_num(k,m) = ISJC.num(ij);
        ISJC2.pos_Chr(k,m) = ISJC.pos2(ij);
        ISJC2.pos_IS(k,m) = ISJC.pos1(ij);
        ISJC2.dir_Chr(k,m) = ISJC.dir2(ij);
        ISJC2.dir_IS(k,m) = ISJC.dir1(ij);
        ISJC2.refIS(k,m) = ISJC.refIS(ij);
        ISJC2.read_count_offset(k,m) = ISJC.read_count_offset(ij);
        ISJC2.IS_distance{k,1} = ISJC.IS_dis{i}(ii)';
        ISJC2.IS_distance{k,2} = ISJC.IS_dis{j}(jj)';
        ISJC2.breseq_reject(k,m) = breseq_reject(ij);

        % Store shared properties of the two junctions:
        ISJC2.scaf_Chr(k) = ISJC.jscaf2(i);
        ISJC2.IS{k} = arrayfun(@directed_IS, IS_id, orientation);
        ISJC2.span_origin(k) = span_origin;
    end
end
end

