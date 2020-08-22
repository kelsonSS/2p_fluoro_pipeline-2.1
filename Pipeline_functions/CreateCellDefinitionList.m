function needs_CellIDs = CreateCellDefinitionList(savepath)

if ~exist('savepath','var')
    savepath = '\\vault3\data\kelson\analyzed';
end 


 registered = dir( fullfile(savepath, '\**\GreenChannelRegistered*.raw'));
 registered = struct2cell(registered);
 registered = registered(2,:)';
 
 CellDefinitions = dir( [savepath, '\**\CellDefinitions*.mat']);
 CellDefinitions = struct2cell(CellDefinitions);
 CellDefinitions = CellDefinitions(2,:)';
 
 
 needs_CellIDs =  setdiff(registered,CellDefinitions);


end 
