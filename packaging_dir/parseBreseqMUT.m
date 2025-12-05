function parseBreseqMUT(infile,outpth)
%read breseq output (output.gd)

if exist(outpth,'file')
    return
end

% Build empty tables
table_names = {'SNP', 'MOB', 'DEL', 'JC','UN'};

fields = struct;
tables = struct;
index = struct;
for i = 1:numel(table_names)
    nm = table_names{i};
    tbl = readtable([get_folder('fields') filesep nm '_fields.csv']);
    fields.(nm) = tbl;
    tables.(nm) = table('Size',[1e4,height(tbl)], ...
        'VariableTypes',tbl.types,'VariableNames',tbl.fields);
    index.(nm) = 0;
end

% Parse BRESEQ into tables:
fin = fopen(infile,'r') ;
while ~feof(fin)
    l = fgetl(fin);
    category = l(1:find([l, char(9)]==char(9),1)-1);
    if ~ismember(category,table_names)
        continue
    end
    tbl = tables.(category);
    flds = fields.(category);
    ind = index.(category)+1;
    n_non_optional = find(~flds.optional ,1, 'last');
    tabs = [0, strfind(l,char(9)), length(l)+1];
    for j = 1:numel(tabs)-1
        str = l(tabs(j)+1:tabs(j+1)-1);
        if j<=n_non_optional
            fld = flds.fields{j};
        else
            parts = split(str, '=');
            fld = parts{1};
            str = parts{2};
        end
        [~,k] = ismember(fld, flds.fields);
        if k==0
            logger(sprintf('Field "%s" appear in BRESEQ output.gd and not listed in Matlab tables.', fld),'ERROR')
        end
        value = convert_str_to_type(str, flds.types{k});
        if isempty(value) 
            if strcmp(str,'NA')
                value = nan ; 
            else
                error('could not parse value for category %s')
            end
        end
            
        tbl.(fld)(ind) = value;
    end
    index.(category) =  index.(category)+1;
    tables.(category) = tbl;
end
fclose (fin) ;

% Remove empty reserve lines from each table:
for i = 1:numel(table_names)
    nm = table_names{i};
    k = index.(nm);
    tables.(nm) = tables.(nm)(1:k,:);
end

% Save:
save(outpth, '-struct', 'tables')

end


function value = convert_str_to_type(str, type)
switch type
    case 'string'
        value = str;
    case 'cell'
        value = {str2num(str)};
    case 'double'
        value = str2num(str);
    case 'int32'
        value = int32(str2num(str));
end
end
