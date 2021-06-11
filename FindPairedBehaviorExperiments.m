function [paired_expts_list,unpaired] = FindPairedBehaviorExperiments()

paired_expts_list = {};
unpaired = {};
active_files = dir('Z:\Kelson\Analyzed\**\*\Tones*Active\');
folders = unique({active_files.folder}');
%to_click = {}
expt_idx = 1; 
for fl_idx = 1:length(folders)
    
    %  find paired passive file if it exists 
    passive_dir = dir( fullfile(...
        fileparts(folders{fl_idx}),...
                  'Tones*Passive',...
                   'Fluorescence.mat'));
    
    if ~isempty(passive_dir)
        passive_file =  FullFileNameFromDir(passive_dir);
    else 
        unpaired{end+1} = folders{fl_idx};
    end 
    active_file = folders{fl_idx};
    
   % add to experiment file and increment idx
    
    paired_expts_list{expt_idx,1} = fileparts(passive_file);
    paired_expts_list{expt_idx,2} = fileparts(active_file);
    
   expt_idx = expt_idx +1 ; 
end 

    
end 

    
    
  
  
  