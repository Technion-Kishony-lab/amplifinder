function folder = get_folder(name)

global baseDir

switch name
    case 'genomesDB'
        folder = 'genomesDB';
    case 'fasta'
        folder = ['genomesDB' filesep 'fasta'];
    case 'ISfinder'
        folder = ['genomesDB' filesep 'ISfinder'];
    case 'genbank'
        folder = ['genomesDB' filesep 'genbank'];
    case 'prms'
        folder = 'prms' ;
    case 'master_prms'
        folder = baseDir ;
    case 'fields'
        folder = fullfile(fileparts(baseDir), 'amplifinder', 'data', 'fields');
end

