function paired_expts_list,to_click = FindPairedBehaviorExperiments()

paired_expts_list = {};

files = dir('Z:\Kelson\Analyzed\**\2020*\Tones*Active\');
to_click = {}
expt_idx = 1; 
for fl_idx = 1:length(files)
    
    %  find paired passive file if it exists 
    passive_dir = dir( fullfile(...
        fileparts(files(fl_idx).folder),...
                  'Tones*Passive',...
                   'Fluorescence.mat'));
    
    if ~isempty(passive_dir)
        passive_file =  FullFileNameFromDir(passive_dir);
    else 
        continue
    end 
    active_file = FullFileNameFromDir(dir(fullfile(files(fl_idx),'Fluorescence.mat');
    
   % add to experiment file and increment idx
    
    paired_expts_list{expt_idx,1} = fileparts(passive_file);
    paired_expts_list{expt_idx,2} = fileparts(active_file);
    
   expt_idx = expt_idx +1 ; 
end 

    
end 

    
    
  
  
  