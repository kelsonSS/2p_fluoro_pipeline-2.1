function  FRA_list = findAnimalFRAS(TNFullBehaviorTable)
analyzed_basefolder = 'Z:\Kelson\Analyzed'

unique_names = unique(TNFullBehaviorTable.AnimalID);

FRA_list = {}

for file_idx = 1:length(unique_names)
    files_found = dir(fullfile(analyzed_basefolder,unique_names{file_idx},'**','FRA*oise'));
    
    folder_names = {};
    if ~isempty(files_found)
        file_names = {files_found.name};
        file_folders ={files_found.folder};
        for folder_idx = 1:length(file_folders)
            folder_names(folder_idx) = fullfile(file_folders(folder_idx),file_names{folder_idx});
            
        end
    end 
    
     FRA_list{file_idx} = folder_names;   
    

end

FRA_list = [unique_names, FRA_list'];


 