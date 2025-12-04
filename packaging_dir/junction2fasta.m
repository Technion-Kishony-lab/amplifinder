function junction2fasta(ISJC2,plus_anc,jc_file,jc_vec)

if nargin<4
    jc_vec = 1:7;
end

for i = 1:height(ISJC2)
    isjc2 = ISJC2(i,:);
    iso_outpath = isjc2.iso_outpath;
    anc_outpath = isjc2.anc_outpath;
    dn = isjc2.analysis_directory;
    if exist(fullfile(iso_outpath,dn,jc_file),'file')
        % for the unusual case of the isolate having the directory and the
        % ancestor not having it
        if plus_anc && ~exist(fullfile(anc_outpath,dn,jc_file),'file')
            mkdir(fullfile(anc_outpath,dn))
            copyfile(fullfile(iso_outpath, dn, jc_file),fullfile(anc_outpath, dn))
        end
        continue
    end

    warning ('off','bioinfo:fastawrite:AppendToFile')

    load (iso_outpath+filesep+'fastq_read_length','read_length')
    load (iso_outpath+filesep+'IS_loc','IS_loc')
    load (iso_outpath+filesep+'ref_props','ref_props')

    scaf_Chr = isjc2.scaf_Chr ;
    Chr_seq = fastaread([char(isjc2.ref_path) '/fasta/' ref_props.name{scaf_Chr} '.fasta']) ;
    Chr_seq = Chr_seq.Sequence ;

    f = isjc2.chosen_IS{1}.id; % the line number of the chosen IS in IS_loc.
    jscaf = IS_loc.IS_jscaf(f);
    IS_seq = fastaread([char(isjc2.ref_path) '/fasta/' ref_props.name{jscaf} '.fasta']);
    IS_seq = IS_seq.Sequence;

    WID = ceil(mean(read_length))*2 ;

    if exist(fullfile(anc_outpath,dn,jc_file),'file') && plus_anc
        if ~exist(fullfile(iso_outpath,dn,jc_file),'file')
            mkdir(fullfile(iso_outpath,dn))
            copyfile(fullfile(anc_outpath,dn,jc_file),fullfile(iso_outpath, dn))
        end
        continue
    end

    mkdir(fullfile(iso_outpath,dn))

    % adapt to origin spanning
    if ~isjc2.span_origin
        ij = [1 2] ;
    else
        ij = [2 1] ;
    end
    pos = isjc2.pos_Chr(ij);
    refIS = isjc2.refIS(ij);
    paired_single_locus_pos = isjc2.pos_of_paired_single_locus_junction(ij);
    IS_orientation = isjc2.chosen_IS{1}.orientation;
    IS_distance(1) = isjc2.IS_distance{ij(1)}(isjc2.IS{1}==isjc2.chosen_IS{1});
    IS_distance(2) = isjc2.IS_distance{ij(2)}(isjc2.IS{1}==isjc2.chosen_IS{1});
    paired_single_locus_dis = zeros(1,2);


    if ~isempty(isjc2.IS_of_paired_single_locus_junction{ij(1)}) && ismember(isjc2.chosen_IS{1},isjc2.IS_of_paired_single_locus_junction{ij(1)})
        paired_single_locus_dis(1) = isjc2.dis_of_paired_single_locus_junction{ij(1)}(isjc2.IS_of_paired_single_locus_junction{ij(1)}==isjc2.chosen_IS{1});
    else
        paired_single_locus_dis(1) = 0;
    end


    if ~isempty(isjc2.IS_of_paired_single_locus_junction{ij(2)}) && ismember(isjc2.chosen_IS{1},isjc2.IS_of_paired_single_locus_junction{ij(2)})
        paired_single_locus_dis(2) = isjc2.dis_of_paired_single_locus_junction{ij(2)}(isjc2.IS_of_paired_single_locus_junction{ij(2)}==isjc2.chosen_IS{1});
    else
        paired_single_locus_dis(1) = 0;
    end

    % adjusting chromosome positions outside of amplicon
    if paired_single_locus_pos(1) > 0
        if refIS(1) > 0
            posOut(1) = pos(1) - (IS_loc.LocRight(f) - IS_loc.LocLeft(f) + 1 ) - 1 ;
            posOutwithIS(1) = pos(1) - (IS_loc.LocRight(f) - IS_loc.LocLeft(f) + 1 ) - 1;
        else
            posOut(1) =  pos(1) - 1 ;
            posOutwithIS(1) = paired_single_locus_pos(1) ;
        end
    else
        posOut(1) = pos(1) - 1 ;
        posOutwithIS(1) = pos(1) - 1 ;
    end
    if paired_single_locus_pos(2) > 0
        if refIS(2) > 0
            posOut(2) = pos(2) + (IS_loc.LocRight(f) - IS_loc.LocLeft(f) + 1 ) + 1;
            posOutwithIS(2) = pos(2) + (IS_loc.LocRight(f) - IS_loc.LocLeft(f) + 1 ) + 1;
        else
            posOut(2) = pos(2) + 1 ;
            posOutwithIS(2) = paired_single_locus_pos(2) ;
        end
    else
        posOut(2) = pos(2) + 1 ;
        posOutwithIS(2) = pos(2) + 1 ;
    end

    % adjusting IS positions

    ISLeft = IS_loc.LocLeft(f) + IS_distance(2);
    ISRight = IS_loc.LocRight(f) - IS_distance(1);

    ISLeft_edge = IS_loc.LocLeft(f) + paired_single_locus_dis(2);
    ISRight_edge = IS_loc.LocRight(f) - paired_single_locus_dis(1);

    %overall view:
    %    ~~~>>>======>>>======>>>~~~
    %(1)  ~~==
    %(2)  ~~>>
    %(3)           ==>>
    %(4)            ====
    %(5)             >>==
    %(6)                       >>~~
    %(7)                       ==~~

    %(1) left reference
    %~~==
    seq1 = Chr_seq(posOut(1)-WID+1:posOut(1)) ;
    seq2 = Chr_seq(pos(1):pos(1)+WID-1) ;
    seq = [seq1,seq2] ;
    if ismember(1,jc_vec)
        fastawrite(fullfile(iso_outpath,dn,jc_file),'1',seq)
    end

    %(2) left IS (trans)
    %~~>>
    seq1 = Chr_seq(posOutwithIS(1)-WID+1:posOutwithIS(1)) ;
    if IS_orientation==1
        seq2 = IS_seq(ISLeft_edge:ISLeft_edge+WID-1) ;
    else
        seq2 = seqrcomplement(IS_seq(ISRight_edge-WID+1:ISRight_edge)) ;
    end
    seq = [seq1,seq2];
    if ismember(2,jc_vec)
        fastawrite(fullfile(iso_outpath,dn,jc_file),'2',seq)
    end

    %(3) left of mid IS
    %==>>
    seq1 = Chr_seq(pos(2)-WID+1:pos(2)) ;
    if IS_orientation==1
        seq2 = IS_seq(ISLeft:ISLeft+WID-1) ;
    else
        seq2 = seqrcomplement(IS_seq(ISRight-WID+1:ISRight)) ;
    end
    seq = [seq1,seq2];
    if ismember(3,jc_vec)
        fastawrite(fullfile(iso_outpath,dn,jc_file),'3',seq)
    end

    %(4) lost IS
    %====
    seq1 = Chr_seq(pos(2)-WID+1:pos(2));
    seq2 = Chr_seq(pos(1):pos(1)+WID-1) ;
    seq = [seq1,seq2];
    if ismember(4,jc_vec)
        fastawrite(fullfile(iso_outpath,dn,jc_file),'4',seq)
    end

    %(5) right of mid IS
    %>>==
    if IS_orientation==1
        seq1 = IS_seq(ISRight-WID+1:ISRight);
    else
        seq1 = seqrcomplement(IS_seq(ISLeft:ISLeft+WID-1)) ;
    end
    seq2 = Chr_seq(pos(1):pos(1)+WID-1);
    seq = [seq1,seq2];
    if ismember(5,jc_vec)
        fastawrite(fullfile(iso_outpath,dn,jc_file),'5',seq)
    end

    %(6) right IS (trans)
    %>>~~
    if IS_orientation==1
        seq1 = IS_seq(ISRight_edge-WID+1:ISRight_edge);
    else
        seq1 = seqrcomplement(IS_seq(ISLeft_edge:ISLeft_edge+WID-1));
    end
    seq2 = Chr_seq(posOutwithIS(2):posOutwithIS(2)+WID-1);
    seq = [seq1,seq2];
    if ismember(6,jc_vec)
        fastawrite(fullfile(iso_outpath,dn,jc_file),'6',seq)
    end

    %(7) right reference
    %==~~
    seq1 = Chr_seq(pos(2)-WID+1:pos(2));
    seq2 = Chr_seq(posOut(2):posOut(2)+WID-1) ;
    seq = [seq1,seq2];
    if ismember(7,jc_vec)
        fastawrite(fullfile(iso_outpath,dn,jc_file),'7',seq)
    end

    if plus_anc
        mkdir(fullfile(anc_outpath,dn))
        copyfile(fullfile(iso_outpath, dn, jc_file),fullfile(anc_outpath, dn))
    end

end