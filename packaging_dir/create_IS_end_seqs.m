function create_IS_end_seqs(ref_name,ref_path,in_span, out_span)
% create a file containing the

source = {'genbank', 'ISfinder'};
for k = 1:2
    outfile = fullfile(char(ref_path),source{k},[ref_name,'_IS_seqs.mat']);
    if exist(outfile, 'file')
        continue
    end
    if ~exist('seq','var')
        fasta = fastaread(fullfile(char(ref_path),'fasta',[ref_name,'.fasta']));
        seq = fasta.Sequence;
    end

    clear IS_loc
    infile = fullfile(char(ref_path),source{k},[ref_name,'_IS_loc.mat']);
    load(infile, 'IS_loc')

    assert(all(IS_loc.IS_scaf == ref_name))
    [IS_seqs, IS_seqs_with_margins, IS_end_seqs] = IS_loc_to_seqs(IS_loc, seq);
    save(outfile, 'IS_seqs','IS_seqs_with_margins','IS_end_seqs','in_span','out_span')
end


    function [IS_seqs,IS_seqs_with_margins,IS_end_seqs] = IS_loc_to_seqs(IS_loc, seq)

    h = height(IS_loc);
    IS_end_seqs = cell(h, 2);
    IS_seqs_with_margins = cell(h, 2);
    IS_seqs = cell(h, 1);
    % to cover the option of a non-circular scaffold ending or starting
    % with an IS element not leaving room for spanning regions. 
    % a fake false region is generated (polyA). 
    elongation_dist = max(in_span,out_span);
    seq = [repmat('a',1,elongation_dist),seq,repmat('a',1,elongation_dist)];
    IS_loc.LocLeft = IS_loc.LocLeft+elongation_dist;
    IS_loc.LocRight = IS_loc.LocRight+elongation_dist;
    for i = 1:h
        l = IS_loc.LocLeft(i);
        r = IS_loc.LocRight(i);
        IS_seqs{i,1} = seq(l:r);
        IS_seqs_with_margins{i,1} = seq(l-out_span:r+out_span) ; 
        IS_seqs_with_margins{i,2} = seqrcomplement(seq(l-out_span:r+out_span)) ; 
        IS_end_seqs{i,1} = seq(l-out_span:l+in_span);
        IS_end_seqs{i,2} = seqrcomplement(seq(r-in_span:r+out_span));
    end
    end

end