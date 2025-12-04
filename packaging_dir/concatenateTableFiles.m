function ConcatenatedTable = concatenateTableFiles(FileList)
ConcatenatedTable = table;
for iFile = 1:numel(FileList)
    Data  = load([FileList{iFile}]);
    Field = struct2cell(Data);
    ConcatenatedTable = [ConcatenatedTable; Field{1}];
end
