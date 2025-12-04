function parse_reads_from_bam(bampath,bamname,outputname)

kd = 0 ;
ki = 0 ;


b = BioMap(fullfile(bampath,bamname)) ;
nmbr_reads = length(b.Signature);
SIGNATURE = b.Signature;
RdStart = [];
RdFlag = [];
RdLength = [];
RdRef = [];
RdFull = [];

DelReg_read = uint32(zeros(nmbr_reads,1));
DelReg_ref = cell(nmbr_reads,1);
DelReg_rdpos = uint32(zeros(nmbr_reads,1));
DelReg_pos = uint32(zeros(nmbr_reads,1));
DelReg_length = uint32(zeros(nmbr_reads,1));
DelReg_length_sup_left = uint32(zeros(nmbr_reads,1));
DelReg_length_sup_right = uint32(zeros(nmbr_reads,1));
DelReg_single = false(nmbr_reads,1);

InsReg_read = uint32(zeros(nmbr_reads,1));
InsReg_ref = cell(nmbr_reads,1);
InsReg_rdpos = uint32(zeros(nmbr_reads,1));
InsReg_pos = uint32(zeros(nmbr_reads,1));
InsReg_length = uint32(zeros(nmbr_reads,1));
InsReg_length_sup_left = uint32(zeros(nmbr_reads,1));
InsReg_length_sup_right = uint32(zeros(nmbr_reads,1));
InsReg_single = false(nmbr_reads,1);

if ~isempty(b.Reference)
    %start
    RdStart = double(b.Start) ;

    %flag
    RdFlag = getFlag(b) ;

    %length
    RdLength = zeros(nmbr_reads,1) ;
    RdFull = false(nmbr_reads,1) ;
    tic
    %referene (1:6)
    RdRef = getReference(b);
    if ~strcmp(class(RdRef{1}),'char')
        RdRef = cellfun(@str2num,RdRef);
    end
    for nn = 1:length(SIGNATURE)
        if ~mod(nn,10000)
            fprintf('%d/%d reads in %.1f\n',nn,nmbr_reads,toc);
            tic
        end
        tst = SIGNATURE{nn} ;
        if length(unique([tst '0' '1' '2' '3' '4' '5' '6' '7' '8' '9' 'M'])) == 11 % e.g. no 'N','S','H','P','I','D'
            tstf1 = strfind(tst,{'M'}) ;
            n = str2num(tst(1:tstf1-1)) ;
            if n==0
                error('length of read == 0')
            else
                RdLength(nn) = n ;
            end
            RdFull(nn) = true ;
        elseif length(unique([tst '0' '1' '2' '3' '4' '5' '6' '7' '8' '9' 'M' 'D' 'I'])) == 13 % e.g. no 'N','S','H','P'
            tstf1 = strfind(tst,{'M'}) ;
            tstf2 = strfind(tst,{'D'}) ;
            tstf3 = strfind(tst,{'I'})  ;
            op = sort([tstf1,tstf2,tstf3]) ;
            if ismember(op(1),tstf1) && ismember(op(end),tstf1)
                op = [0;op'] ;
                for mm = 2:length(op)
                    n = str2num(tst(op(mm-1)+1:op(mm)-1)) ;
                    if strcmp(tst(op(mm)),'M') 
                        RdLength(nn) = RdLength(nn) + n ;
                    end
                    if strcmp(tst(op(mm)),'D')
                        kd = kd + 1 ;
                        DelReg_read(kd) = nn ; % read id
                        DelReg_ref(kd) = RdRef(nn);
                        DelReg_rdpos(kd) = RdLength(nn) ; % position of deletion within the read
                        DelReg_pos(kd) = RdStart(nn) + RdLength(nn) - 1 ; % position of deletion in reference
                        DelReg_length(kd) = n ; % deletion length
                        if mm>2
                            if strcmp(tst(op(mm-1)),'M')
                                DelReg_length_sup_left(kd) = str2num(tst(op(mm-2)+1:op(mm-1)-1)) ; %support for alignment from left
                            else
                                DelReg_length_sup_left(kd) = 0 ;
                            end
                        else
                            DelReg_length_sup_left(kd) = 0 ;
                        end
                        if length(op) > mm
                            if strcmp(tst(op(mm+1)),'M')
                                DelReg_length_sup_right(kd) = str2num(tst(op(mm)+1:op(mm+1)-1)) ; % support for alignment from right
                            else
                                DelReg_length_sup_right(kd) = 0 ;
                            end
                        else
                            DelReg_length_sup_right(kd) = 0 ;
                        end
                        if length(tstf1)==2 && length(op)==4
                            DelReg_single(kd) = true ;
                        else
                            DelReg_single(kd) = false ;
                        end
                        RdLength(nn) = RdLength(nn) + n ;
                    end
                    if strcmp(tst(op(mm)),'I')
                        ki = ki + 1 ;
                        InsReg_read(ki) = nn ;
                        InsReg_ref(ki) = RdRef(nn);
                        InsReg_rdpos(ki) = RdLength(nn) ;
                        InsReg_pos(ki) = RdStart(nn) + RdLength(nn) -1 ;
                        InsReg_length(ki) = n ;
                        if mm>2
                            if strcmp(tst(op(mm-1)),'M')
                                InsReg_length_sup_left(ki) = str2num(tst(op(mm-2)+1:op(mm-1)-1)) ;
                            else
                                InsReg_length_sup_left(ki) = 0 ;
                            end
                        else
                            InsReg_length_sup_left(ki) = 0 ;
                        end
                        if length(op) > mm
                            if strcmp(tst(op(mm+1)),'M')
                                InsReg_length_sup_right(ki) = str2num(tst(op(mm)+1:op(mm+1)-1)) ;
                            else
                                InsReg_length_sup_right(ki) = 0 ;
                            end
                        else
                            InsReg_length_sup_right(ki) = 0 ;
                        end
                        if length(tstf1)==2 && length(op)==4
                            InsReg_single(ki) = true ;
                        else
                            InsReg_single(ki) = false ;
                        end
                    end
               end
           end
        end
    end
    
end
DelReg = table(DelReg_read, DelReg_ref, DelReg_rdpos, DelReg_pos, DelReg_length, DelReg_length_sup_left, DelReg_length_sup_right, DelReg_single,'VariableNames',{'read','ref','rdpos','pos','length','length_sup_left','length_sup_right','single'}) ;
InsReg = table(InsReg_read, InsReg_ref, InsReg_rdpos, InsReg_pos, InsReg_length, InsReg_length_sup_left, InsReg_length_sup_right, InsReg_single,'VariableNames',{'read','ref','rdpos','pos','length','length_sup_left','length_sup_right','single'}) ;
save (strcat(bampath, filesep, outputname),'RdStart','RdFlag','RdLength','RdRef','RdFull','DelReg','InsReg')