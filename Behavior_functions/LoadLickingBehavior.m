function Behavior = LoadLickingBehavior(Folders)


if iscell(Folders)
    n_folders = length(Folders);
else 
    n_folders = 1;
end 
    
Behavior = cell(n_folders,1);
for file_idx = 1:n_folders

    lick_file_path = fullfile(Folders{file_idx},'BehavioralResponses_Frames.mat');
    Behavior(file_idx) = {tryLoading(lick_file_path)};
end 

end


function LickTimes  = tryLoading(file)

LickTimes= {};
try 
LickTimes = load(file);
catch 
end 

end 
