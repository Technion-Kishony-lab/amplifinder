function ISJC2 = calculate_amplicon_length(ISJC2, ref_props)
% Calculate the distance between junction pairs, accounting for their
% direction and circularity of the chromosme. 

global prms

% Our junction pairs are facing opposite directions on the chromosome:
assert(all(prod(ISJC2.dir_Chr,2) == -1));

is_circular = ref_props.circular(ISJC2.scaf_Chr);
scaf_length = int32(ref_props.length(ISJC2.scaf_Chr));
span_origin = ISJC2.span_origin; 
complementary_length = zeros(height(ISJC2),1);

% -----|=======|----
amplicon_length = diff(ISJC2.pos_Chr, [], 2)+1;


% Circular chromosome:
% =====|-------|====
%   <--         -->
z1 = span_origin & is_circular;
z2 = ~span_origin & is_circular; 
complementary_length(z1) = amplicon_length(z1);
complementary_length(z2) = scaf_length(z2) - amplicon_length(z2);
amplicon_length(z1) = scaf_length(z1) - amplicon_length(z1);

% Linear chromosome:
% -----|-------|----
%   <--         -->
z1 = span_origin & ~is_circular;
z2 = ~span_origin & ~is_circular;
complementary_length(z1) = amplicon_length(z1);
complementary_length(z2) = inf;
amplicon_length(z1) = inf;

ISJC2.amplicon_length = amplicon_length;
ISJC2.complementary_length = complementary_length; 
end