function [RelativePath] = RelativePathFromFullPath(FullPath,in_path)
%


  n = length(in_path);
  
  % remove base path
  RelativePath = FullPath(n+1:end);