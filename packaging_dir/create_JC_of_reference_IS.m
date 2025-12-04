function JC_refIS = create_JC_of_reference_IS(iso_outpath)
% parse ancestral IS elements into JC format

global prms
REFERENCE_IS_OUT_SPAN = prms.REFERENCE_IS_OUT_SPAN;  % extracting a sequence, long enough to be unique in the chromosome (100)

filename = fullfile(iso_outpath, 'JC_refIS.mat');

if exist(filename, 'file')
    load(filename, 'JC_refIS')
    return
end

load(fullfile(iso_outpath, 'IS_loc.mat'), 'IS_loc')

JC_refIS = table;
for i = 1:height(IS_loc)
    % for each IS we make two lines representing the two ends of 
    % IS-chromosomal junctions:
    reference_IS_in_span = IS_loc.LocRight(i)-IS_loc.LocLeft(i)+1; 
    tmp.type(1:2,1) = ["JC";"JC"] ;
    tmp.num(1:2,1) = [0;0] ;
    tmp.dot(1:2,1) = [".";"."] ;
    tmp.jscaf1(1:2,1) = [IS_loc.IS_jscaf(i);IS_loc.IS_jscaf(i)] ;
    tmp.pos1(1:2,1) = int32([IS_loc.LocLeft(i);IS_loc.LocRight(i)]) ;
    tmp.dir1(1:2,1) = [1;-1] ;
    tmp.jscaf2(1:2,1) = [IS_loc.IS_jscaf(i);IS_loc.IS_jscaf(i)] ;
    tmp.pos2(1:2,1) = int32([IS_loc.LocLeft(i) - 1 ;IS_loc.LocRight(i) + 1]) ; 
    tmp.dir2(1:2,1) = [-1;1] ;
    tmp.refIS(1:2,1) = [i;i] ; 
    tmp.zero(1:2,1) = [0;0] ;
    tmp.alignment_overlap(1:2,1) = [0;0] ;
    tmp.circular_chromosome(1:2,1) = [0;0] ;
    tmp.coverage_minus(1:2,1) = [0;0] ;
    tmp.coverage_plus(1:2,1) = [0;0] ;
    tmp.flanking_left(1:2,1) = int32([reference_IS_in_span;REFERENCE_IS_OUT_SPAN]) ;   %**
    tmp.flanking_right(1:2,1) = int32([reference_IS_in_span;REFERENCE_IS_OUT_SPAN]) ;  %**
    tmp.frequency(1:2,1) = [0;0] ;
    tmp.junction_possible_overlap_registers(1:2,1) = [0;0] ;
    tmp.key(1:2,1) = ["";""] ;
    tmp.max_left(1:2,1) = [0;0] ;
    tmp.max_left_minus(1:2,1) = [0;0] ;
    tmp.max_left_plus(1:2,1) = [0;0] ;
    tmp.max_min_left(1:2,1) = [0;0] ;
    tmp.max_min_left_minus(1:2,1) = [0;0] ;
    tmp.max_min_left_plus(1:2,1) = [0;0] ;
    tmp.max_min_right(1:2,1) = [0;0] ;
    tmp.max_min_right_minus(1:2,1) = [0;0] ;
    tmp.max_min_right_plus(1:2,1) = [0;0] ;
    tmp.max_pos_hash_score(1:2,1) = [0;0] ;
    tmp.max_right(1:2,1) = [0;0] ;
    tmp.max_right_minus(1:2,1) = [0;0] ;
    tmp.max_right_plus(1:2,1) = [0;0] ;
    tmp.neg_log10_pos_hash_p_value(1:2,1) = [0;0] ;
    tmp.new_junction_coverage(1:2,1) = [0;0] ;
    tmp.new_junction_read_count(1:2,1) = [0;0] ;
    tmp.polymorphism_frequency(1:2,1) = [0;0] ;
    tmp.pos_hash_score(1:2,1) = [0;0] ;
    tmp.prediction(1:2,1) = ["";""] ;
    tmp.reject(1:2,1) = ["";""] ;
    tmp.side_1_annotate_key(1:2,1) = ["";""] ;
    tmp.side_1_continuation(1:2,1) = [0;0] ;
    tmp.side_1_coverage(1:2,1) = [0;0] ;
    tmp.side_1_overlap(1:2,1) = [0;0] ;
    tmp.side_1_possible_overlap_registers(1:2,1) = [0;0] ;
    tmp.side_1_read_count(1:2,1) = [0;0] ;
    tmp.side_1_redundant(1:2,1) = [0;0] ;
    tmp.side_2_annotate_key(1:2,1) = ["";""] ;
    tmp.side_2_continuation(1:2,1) = [0;0] ;
    tmp.side_2_coverage(1:2,1) = [0;0] ;
    tmp.side_2_overlap(1:2,1) = [0;0] ;
    tmp.side_2_possible_overlap_registers(1:2,1) = [0;0] ;
    tmp.side_2_read_count(1:2,1) = [0;0] ;
    tmp.side_2_redundant(1:2,1) = [0;0] ;
    tmp.total_non_overlap_reads(1:2,1) = [0;0] ;
    tmp.unique_read_sequence(1:2,1) = ["";""] ;
    tmp.read_count_offset(1:2,1) = [0;0] ;

    JC_refIS = [JC_refIS;struct2table(tmp)];
end

save(filename, 'JC_refIS')
