function [Active_list,experiment_files] = Fluoro_to_Table_Behavior(experiment_files)
% This function extends Fluoro_to_Table to work with behavioral experiments
%
% since behavioral experiments can run for a variable number of trials
% depending on the animal, we cannot concatenate files from different
% experiments together. the solution is to iteratively call Fluoro_to_table
% for each file and then put them together in a cell

% in this way you can interact with active files the same as passive files
% with the addition of the bracket indicating the experiment you want to
% sample from 
% ex. active{1}.DFF = the DFF from the first active experiment  


% Kelson Shilling-Scrivo -2019 
if ~exist('experiment_files','var')
    run = 1;
    expt = 1;
    while run 
        experiment_files{expt,1} = uigetdir('\\vault3\data\kelson\Analyzed','Load Passive')
        experiment_files{expt,2} = uigetdir(experiment_files{expt,1},'Load Active')
        expt = expt+1;
        
       build_experiments_flag =  questdlg('continue?');
       if ~strcmp(build_experiments_flag, 'Yes')
           run = 0 ;
       end 
    end 
end 

    

m = length(experiment_files);
 
 
Active_list=  cell(length(experiment_files),1) ; 

for expt_id = 1:m
    try
        % get passive data
    passive = Fluoro_to_Table(experiment_files{expt_id,1});
    Active_list{expt_id,1} = passive;
    Active_list{expt_id,3} = experiment_files{expt_id,1};
    
  
    catch
 
    end 
    try
        % get active data
    active =  Fluoro_to_Table(experiment_files{expt_id,2});  
    Active_list{expt_id,2} = active; 
    Active_list{expt_id,4} = experiment_files{expt_id,2};
    catch
        
    end 
    
end 

