function[paths,expname,psignalfiles, matless,animalIDs] = createDataList(file_path,all_files)

% this function checks the folder containing every raw file in 
% file_path and checks to see if there is a corresponding 
% folder in file_path\Analyzed. if there is not that indicates that the
% file needs to be analyed and creates the appropriate Datalist.
%
% all_files - 0 to return data which has not yet finished the pipeline 
%              (has a fluorescence.mat file)
%           -  1 to return all data in the folder  
%
% This file is part of the 2p_fluoro_pipeline 2.1


if ~exist('file_path','var')
file_path = '\\VAULT3\Data\Kelson\Files to unload' ;
end 

analyzed_path = fullfile(file_path,'analyzed');

 unregistered = dir( [file_path, '\**\Image0001_0001.raw']);
 unregistered = struct2cell(unregistered);
 unregistered = unregistered(2,:)';

 
analyzed_idx = cellfun(@isanalyzed,unregistered,0)
 
% subset to only incomplete files 
if ~all_files
    unregistered = unregistered(~analyzed_idx);
end 







    function analyzed_list = checkIfAnalyzed(files_to_check,in_path,analyzed_path)
        
        analyzed_path_to_check = strrep(files_to_check,in_path,analyzed_path)
        
        analyzed_fluoro_path = fullfile(analyzed_path_to_check,...
                                        'Fluorescence.mat')
        
        analyzed_list = exist('analyzed_Fluoro_path','file');
    end 
        
        
        
        
        
        
        
        
        
        
        
 
 


end     
    
 
     

