function needs_Extraction = CreateCellExtractionList(savepath)

if ~exist('savepath','var')
    savepath = '\\vault3\data\kelson\analyzed';
end 

 
 CellDefinitions = dir( [savepath, '\**\CellDefinitions.mat']);
 CellDefinitions = struct2cell(CellDefinitions);
 CellDefinitions = CellDefinitions(2,:)';
 
 Fluorescence = dir( [savepath, '\**\Fluorescence.mat']);
 Fluorescence = struct2cell(Fluorescence);
 Fluorescence = Fluorescence(2,:)';
 
 
 needs_Extraction =  setdiff(CellDefinitions,Fluorescence);


end 
