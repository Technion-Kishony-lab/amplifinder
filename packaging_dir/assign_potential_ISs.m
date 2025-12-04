function JC = assign_potential_ISs(JC,ref_path,ref_names,isfinder,max_dist_to_IS)
% Find potential IS element ends for each JC

IS_seqs = get_IS_seqs(isfinder, ref_path, ref_names, max_dist_to_IS);

JC_seqs = get_JC_seqs(JC,ref_path,ref_names);

switched_v = false(height(JC),1);
is_ISJC = false(height(JC),1);
for i = 1:height(JC)
    [match1, pos1] = cellfun(@(ISseq)contains_at_start(ISseq, JC_seqs{i,1}, max_dist_to_IS), IS_seqs);
    [match2, pos2] = cellfun(@(ISseq)contains_at_start(ISseq, JC_seqs{i,2}, max_dist_to_IS), IS_seqs);
    is_side1IS = any(match1, [1,2]); 
    is_side2IS = any(match2, [1,2]);
    
    %JCs are 'aligned such that the first 'part' of the junction is the one
    %attributed to an IS
    switched = ~is_side1IS & is_side2IS;
    if switched
        match = match2;
        pos = pos2;
    else
        match = match1;
        pos = pos1;
    end
    is_ISJC(i) = is_side1IS|is_side2IS;
    [f, jside] = find(match);
    ind = sub2ind(size(match), f, jside);
    JC.IS_side{i} = arrayfun(@IS_side, f', jside'*2 - 3);
    JC.IS_dis{i} = int32(pos(ind)' - (max_dist_to_IS+1));
    switched_v(i) = switched;
end

% Switch BRESEQ "1"/"left" <--> "2"/"right", 
% so that "1"/"left" is the side of the IS:
JC.switched = switched_v;
JC(JC.switched,:) = switch_fields(JC(JC.switched,:));

% Remove junctions not involving IS or involving two IS:
JC = JC(is_ISJC, :);

% Sort by scaffold and position:
JC = sortrows(JC,{'jscaf2','pos2'}) ;

end


function [tf, pos] = contains_at_start(str, pat, threshold)
f = strfind(str, pat);
if isempty(f)
    pos = nan;
    tf = false;
else
    pos = f(1);
    tf = f < threshold*2;
end
end

function IS_sequences = get_IS_seqs(isfinder, ref_path, ref_names, max_dist_to_IS)

% Get IS_loc either from genbank or from ISfinder:
folders = {'genbank','ISfinder'};
folder = folders{isfinder+1};

% load IS sequences:
IS_sequences = cell(0, 2); %[left_seq, right_seq]
for i = 1:numel(ref_names)
    load(fullfile(ref_path,folder,[ref_names{i},'_IS_seqs.mat']), ...
        'IS_seqs_with_margins','out_span')
    IS_sequences = [IS_sequences; IS_seqs_with_margins];
    assert(out_span == max_dist_to_IS)
end

end


function seqs = get_JC_seqs(JC,ref_path,ref_names)
% returns the sequences to the two ends of each junction

TRIM_JC_FLANKING_LENGTH = 5 ; %trimming the edges of junction region to avoid misalignments

for i = 1:numel(ref_names)
    % get fasta for all scaffolds of the isolate genome
    fasta = fastaread(fullfile(char(ref_path),'fasta',[ref_names{i} '.fasta']));
    scaf_seqs{i} = fasta.Sequence;
end

seqs = cell(height(JC),2);
for i = 1:height(JC)
    flanking_length = JC.flanking_left(i)-TRIM_JC_FLANKING_LENGTH-1;
    scaf_seq = scaf_seqs{JC.jscaf1(i)};
    seq1 = get_sequence_from_scaffold(scaf_seq, JC.pos1(i), flanking_length, JC.dir1(i));

    flanking_length = JC.flanking_right(i)-TRIM_JC_FLANKING_LENGTH-1;
    scaf_seq = scaf_seqs{JC.jscaf2(i)};
    seq2 = get_sequence_from_scaffold(scaf_seq, JC.pos2(i), flanking_length, JC.dir2(i));

    seqs(i,:) = [{seq1}, {seq2}];
end
end

function seq = get_sequence_from_scaffold(scaf_seq, start, len, dir)
if dir==1
    seq = scaf_seq(start + (0:len));
else
    seq = scaf_seq(start + (-len:0));
    seq = seqrcomplement(seq);
end
end


function JC = switch_fields(JC)
% Swap 1 <--> 2 fields in JC

pairs = {
    'jscaf1'                                'jscaf2'
    'pos1'                                  'pos2'
    'dir1'                                  'dir2'
    'flanking_left'                         'flanking_right'
    'max_left'                              'max_right'
    'max_left_minus'                        'max_right_minus'
    'max_left_plus'                         'max_right_plus'
    'max_min_left'                          'max_min_right'
    'max_min_left_minus'                    'max_min_right_minus'
    'max_min_left_plus'                     'max_min_right_plus'
    'side_1_annotate_key'                   'side_2_annotate_key'
    'side_1_continuation'                   'side_2_continuation'
    'side_1_coverage'                       'side_2_coverage'
    'side_1_overlap'                        'side_2_overlap'
    'side_1_possible_overlap_registers'     'side_2_possible_overlap_registers'
    'side_1_read_count'                     'side_2_read_count'
    'side_1_redundant'                      'side_2_redundant'
    };

fn = fieldnames(JC);
[~,zz] = ismember(pairs,fn);

JC(:, zz(:,[1 2])) = JC(:, zz(:,[2 1]));

end
