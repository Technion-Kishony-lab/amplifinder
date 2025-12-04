function [nmbr_green_reads,nmbr_left_reads,nmbr_right_reads] = count_reads(outputfile,seq_length,read_length,RdStart,RdRef,RdLength,RdFlag,RdFull)

% counts junction supporting reads and non-junction reads (reported as nmbr_XXX_reads)
% saves counts to file together with randomly picked reads for figure, number of reads per panel is determined by yl)
global prms
scafs = {'1','2','3','4','5','6','7'};

if ~exist(outputfile,'file')
    %**
    yl = 20;  % max y limit

    njcr = 10;  % number of jc reads to plot
    req_ovrlp = prms.REQ_OVERLAP ; % minimal bp overlap of jc to verify
    min_bp_in_frame = 10 ;
   
    okz_length = RdLength > read_length*0.9 & RdLength < read_length*1.1 ;
    okz_left = okz_length & RdStart + RdLength > seq_length./2 - RdLength + min_bp_in_frame & (RdStart+RdLength-1) < seq_length./2 + req_ovrlp ;  %filtering for reads within the plotted frame
    okz_right = okz_length & RdStart > seq_length./2 - req_ovrlp & RdStart < seq_length./2 + RdLength - min_bp_in_frame  ;
    okz_green = okz_length & RdStart < seq_length./2 - req_ovrlp & (RdStart+RdLength-1) > seq_length./2 + req_ovrlp ;

    % choosing the reads, green, left, right
    for j = 1:length(scafs)
        ok_green = find(okz_green & RdFull & ismember(RdRef,scafs(j)));
        nmbr_green_reads(1,j) = length(ok_green) ;
        if length(ok_green) > njcr
            temp_rand_perm = sort(randperm(length(ok_green),njcr)) ;
            ok_green = ok_green(temp_rand_perm);
        end
        green_reads.start{j} = RdStart(ok_green) ;
        green_reads.length{j} = RdLength(ok_green) ;
        green_reads.flag{j} = RdFlag(ok_green) ;
        nnjcr = ceil((yl - length(green_reads.start{j}))/2) ; % number of non-jc reads
        ok_left = find(okz_left & RdFull & ismember(RdRef,scafs(j)));
        nmbr_left_reads(1,j) = length(ok_left) ;

        if length(ok_left) > nnjcr
            temp_rand_perm = sort(randperm(length(ok_left),nnjcr)) ;
            ok_left = ok_left(temp_rand_perm);
        end
        left_reads.start{j} = RdStart(ok_left) ;
        left_reads.length{j} = RdLength(ok_left) ;
        left_reads.flag{j} = RdFlag(ok_left) ;

        ok_right = find(okz_right & RdFull & ismember(RdRef,scafs(j))) ;
        nmbr_right_reads(1,j) = length(ok_right) ;
        if length(ok_right) > nnjcr
            temp_rand_perm = sort(randperm(length(ok_right),nnjcr)) ;
            ok_right = ok_right(temp_rand_perm);
        end
        right_reads.start{j} = RdStart(ok_right) ;
        right_reads.length{j} = RdLength(ok_right) ;
        right_reads.flag{j} = RdFlag(ok_right) ;
    end
    save (outputfile, 'green_reads','left_reads','right_reads','nmbr_green_reads','nmbr_left_reads','nmbr_right_reads','yl')
else
    load(outputfile,'nmbr_green_reads','nmbr_left_reads','nmbr_right_reads')
end
