function iso = analyze_fastq_files(iso)

% (1) Generating isolate-specific .mat files in analysis path with fastq stats.
% (2) Read length uniformity is determined. Criterion: sampled reads are all within 5% length difference from each other.
% (3) Adding fastq stats columns to iso:
%     read_length (length of first read in last .fastq file),
%     read_length_uniformity (0/1),
%     num_bases (sum of read_length of each file times number of reads in file).

global prms
sprintf('analyzing read length\n');
for m = 1:height(iso)
    fastq_path = iso.fastq_path(m);

    iso_outpath = char(iso.iso_outpath(m));
    output_file = [iso_outpath filesep 'fastq_size.mat'];

    if exist(output_file,'file')
        load([iso_outpath filesep 'fastq_read_length.mat'],"read_length")
        load([iso_outpath filesep 'fastq_read_length_stats.mat'],"read_length_stats")
        load(output_file,'num_bases')
    else

        % counting total number of reads
        dfgz = dir([fastq_path{1} filesep '*.fastq.gz']) ;
        df = dir([fastq_path{1} filesep '*.fastq']) ;
        assert(isempty(dfgz) | isempty(df))

        if ~isempty(df)
            d = df;
        else
            d = dfgz;
        end

        fldr = {d.folder};
        nm = {d.name};
        fl = strcat(fldr,filesep,nm);
        num_bases = zeros(length(fl),1);
        read_length_stats = cell(1, length(fl));
        read_length = -inf;
        for i = 1:length(fl)
            [num_lines,read_length_temp,read_length_stats1] = count_lines_and_read_length(fl{i});
            num_bases(i) = read_length_temp*num_lines/4;
            read_length_stats{i} = read_length_stats1;
            read_length = max(read_length,read_length_temp);
        end

        save ([iso_outpath filesep 'fastq_read_length.mat'],"read_length") % only read_length of last file is saved!
        save (output_file,"num_bases")
        save ([iso_outpath filesep 'fastq_read_length_stats.mat'],"read_length_stats")
    end

    iso.read_length_uniformity(m) = analyze_read_length_stats(read_length_stats);
    iso.num_bases(m) = sum(num_bases);
    iso.read_length(m) = read_length;
end

if any(isnan(iso.read_length_uniformity))
    warning ('Read length uniformity unknown. Non-uniform length may affect junction coverage accuracy.\n')
elseif any(~iso.read_length_uniformity)
    warning ('Read length is not uniform. May affect junction coverage accuracy.\n')
end

if any(iso.num_bases<prms.MIN_NUM_BASES)
    warning('small .fastq file size');
    %error('.fastq file size is lower than required. Exiting...\n');
end

end

function read_length_uniformity = analyze_read_length_stats(read_length_stats)
tmpv = [read_length_stats{:}];
tmpm = tmpv./tmpv';
if all(tmpm(:)<1.05)   %allowing for a five percent length change between sampled reads
    read_length_uniformity = 1;
elseif any(isnan(tmpm(:)))
    read_length_uniformity = nan;
else
    read_length_uniformity = 0;
end

end
