function needs_CellIDS = CreateCellDefinitionList(savepath)

if ~exist('savepath','var')
    savepath = '\\vault3\data\kelson\analyzed';
end 


 registered = dir( [savepath, '\**\Greenchannelregistered.raw']);
 registered = struct2cell(registered);
 registered = registered(2,:)';
 
 CellDefinitions = dir( [savepath, '\**\CellDefinitions.mat']);
 CellDefinitions = struct2cell(CellDefinitions);
 CellDefinitions = CellDefinitions(2,:)';
 
 
 needs_CellIDs = needs_CellIDS =  setdiff(registered,CellDefinitions);


end 
