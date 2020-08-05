function [RelativePath] = AnimalIDFromPath(FullPath,in_path)
% animal ID should be the first directory after the base path assuming 
% the correct file structure. this code parses the path to obtain the
% animal ID


  n = length(in_path);
  
  % remove base path
  RelativePath = FullPath(n+1:end);