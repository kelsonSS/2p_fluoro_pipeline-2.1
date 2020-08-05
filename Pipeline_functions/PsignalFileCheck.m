function  psignal_files = PsignalFileCheck(Directory) 

psignal_files =[];

PathFiles =  dir(Directory);
PathFiles = struct2cell(PathFiles);
PathFiles = PathFiles(1,:);
% find files
expression = '(ART)|(RND)'; 
matches =  regexp( PathFiles,expression,'once');

% convert logicnl index to row number
psignal_file_idx = find(~cellfun(@isempty,matches));

% extract file 
if psignal_file_idx
    psignal_files = PathFiles{psignal_file_idx};
end 


end 
