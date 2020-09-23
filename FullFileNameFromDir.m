function  fname = FullFileNameFromDir(directory)
% utility function that takes a file result from dir() function and returns
% the full file path 

fname = fullfile(directory.folder,directory.name);

end 

