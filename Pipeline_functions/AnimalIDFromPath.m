function [AnimalID] = AnimalIDFromPath(RelativePath, BasePath)
% animal ID should be the first directory after the base path assuming 
% the correct file structure. this code parses the path to obtain the
% animal ID


if isstr(RelativePath)
    RelativePath = {RelativePath};
end

for ii = 1:length(RelativePath)
    
   CurrPath = RelativePath{ii};

if ~exist('BasePath','var')
    
    if  contains(RelativePath{ii},'Z:\')
        BasePath = 'Z:\Kelson\Analyzed\';
    elseif contains(RelativePath{ii},'\\Vault3')
        BasePath = '\\Vault3\Data\Kelson\Analyzed';
    else
        
        warning('Unknown Base Path. Assuming, Relative Path Given')
        BasePath = '';
    end 
        
end 
    L = length(BasePath);
    CurrPath = CurrPath(L+1:end);
    

parts = strsep(CurrPath,filesep);

parts =  parts(~cellfun(@isempty, parts));

AnimalID{ii} = parts{1};
end 



 