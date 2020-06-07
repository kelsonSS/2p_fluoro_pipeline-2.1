function [experiment_names,animal_names] =  GetExperimentNames(Passive)
%
% helper function to convert old passive experiment files to parse the data
% directory in order to get the proer experiment names 
animal_names ={};
experiment_names = {};
for expt = 1:length(Passive.DataDirs)
    
    
    p = Passive.DataDirs{expt};
    
    [p,expt_type] = fileparts(p);
     [p,date] = fileparts(p);
     [~,animal_name] = fileparts(p);
     
     animal_names{expt} = animal_name;
     
     experiment_names{expt} = [expt_type, '_', animal_name ,'_',date];
    
     
end 
    
    
