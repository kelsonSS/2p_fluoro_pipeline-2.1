function [paired_expts_list,unpaired] = FindPairedBehaviorExperiments()

paired_expts_list = {};
unpaired = {};
active_files = dir('Z:\Kelson\Analyzed\**\*\*ones*ctive\Fluorescence.mat');
folders = unique({active_files.folder}');
%to_click = {}
expt_idx = 1; 
for fl_idx = 1:length(folders)
     
    active_folder = folders{fl_idx};
    
    %  find paired passive file if it exists 
    passive_dir = dir( fullfile(...
        fileparts(active_folder),...
                  'Tones*Passive',...
                   'Fluorescence.mat'));
    
    if ~isempty(passive_dir)
        passive_file =  FullFileNameFromDir(passive_dir);
        passive_folder = fileparts(passive_file);
    else 
        unpaired{end+1} = active_folder;
        continue
    end 
    active_file = folders{fl_idx};
    
   % add to experiment file and increment idx
    
    paired_expts_list{expt_idx,1} = passive_folder;
    paired_expts_list{expt_idx,2} = active_folder;
    
   expt_idx = expt_idx +1 ; 
end 

    
end 

    
    
  
  
  