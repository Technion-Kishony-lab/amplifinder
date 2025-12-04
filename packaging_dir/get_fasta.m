function get_fasta(ref,ref_path)
% downloading fasta from NCBI and simplify header

base_path = [char(ref_path) filesep 'fasta'];
fname = [base_path filesep ref '.fasta']; 
fname_raw = [base_path filesep ref '_raw.fasta'] ;
if ~exist(fname,'file')
    [~,~] = mkdir(char(ref_path));
    [~,~] = mkdir(base_path);
    getgenbank(ref,'ToFile',fname_raw,'FileFormat','FASTA') 
    f = fastaread(fname_raw) ; 
    fastawrite(fname,ref,f.Sequence)
end
