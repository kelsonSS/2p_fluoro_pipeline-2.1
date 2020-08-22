function [AnimalID] = AnimalIDFromPath(RelativePath)
% animal ID should be the first directory after the base path assuming 
% the correct file structure. this code parses the path to obtain the
% animal ID


parts = strsep(RelativePath,filesep);

parts =  parts(~cellfun(@isempty, parts));

AnimalID = parts{1};




 