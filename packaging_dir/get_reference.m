function ref_name = get_reference(input_ref_name,ref_path)
% Create genebank and basic properties of the scaffold input_ref_name:
% creates:
%
% XXXXX - ref_name from NCBI
% RRRRR - input reference name
%
% GENEBANK:
% genomeDB/genebank/RRRRRref.mat: key genome properties
% genomeDB/genebank/XXXXX.mat: matlab genebank
% genomeDB/genebank/XXXXX.GB: NCBI genbank
%
% FASTA:
% genomeDB/fasta/XXXXX.fasta: fasta file with short header (ref_name)
% genomeDB/fasta/XXXXX_raw.fasta: original NCBI fasta
%
% IS:
% genomeDB/genebank/XXXXX_IS_loc_raw.txt: grep of IS lines from XXXXX.GB
% genomeDB/genebank/XXXXX_IS_loc.mat: IS table from XXXXX_IS_loc_raw.txt
%
% genomeDB/ISfinder/XXXXX_blast.txt: ISfinder blast result (IS-db on XXXXX.fasta)
% genomeDB/ISfinder/XXXXX_IS_loc.mat: IS table from blast results
%

ref_name = efetch_genbank(input_ref_name,ref_path);
get_fasta(ref_name,ref_path)
