function[paths,expname,psignalfiles, matless,animalIDs] = createDataList(file_path)

% this function checks the folder containing every raw file in 
% \\VAULT2\Vault2Data\Kelson and checks to see if there is a corresponding 
% folder in ..\Kelson\Analyzed. if there is not that indicates that the
% file needs to be analyed and creates the appropriate Datalist.
%
% This file is part of the 2p_fluoro_pipeline 2.1 and the output variable
% should be called manyDataList 
%
% manyDataList = createDataList
%
% the second optional output matless gives a list of all files that were
% not analyzed because there is not a corresponding behavior (ART) file.
% 
if ~exist('file_path','var')
file_path = '\\VAULT3\Data\Kelson\Files to unload' ;
end 

 % find all unique directories 
 allfiles = dir([file_path , '\**\*.raw']);
 allfiles = struct2cell(allfiles);
 allfiles = allfiles(2,:)';
  % make sure not to count any files that are in analyzed folder
 mask_raw = strfind(allfiles,fullfile(file_path,'Analyzed'));
 mask_raw = cellfun(@isempty,mask_raw);
 file_shorten = @(x) x(length(file_path)+1:end);
            
 raw      = allfiles(mask_raw);
 raw      = cellfun(file_shorten,raw,'UniformOutput',0);
 % find analyzed files by looking for Fluorescence traces
 analyzed = dir([file_path , '\**\Fluorescence.mat']);
 analyzed = struct2cell(analyzed);
 analyzed = analyzed(2,:);
 analyzed = cellfun(file_shorten,analyzed,'UniformOutput',0);
 analyzed = cellfun(@(x) x(10:end),analyzed,'UniformOutput',0);
 
  % also need to show which files need behavior (ART) files so we can 
  % add them if they are missing
  
 matfiles = dir([file_path, '\**\ART*.mat']);
 matfiles = [matfiles;dir([file_path,'\**\RND*.mat'])] ; 
 matfiles = struct2cell(matfiles);
 ARTs     = matfiles(1,:)';
 matfiles = matfiles(2,:)';
 
 matfiles = cellfun(file_shorten,matfiles,'UniformOutput',0);
 % remove all analyzed files
 mask_mat = ~contains(matfiles,'\Analyzed');
 matfiles = matfiles(mask_mat);
 ARTs     = ARTs(mask_mat);       % we only want to look at the raw files 
 
 
 
 
 % test to see if there is an analyzed version for each unanalyzed version 
 % if not add it to the datalist
 unanalyzed = setdiff (raw,analyzed);
 matless    = setdiff(unanalyzed,matfiles); % of those which have behavior files     
 % generate to-do list and associated files 
 todo_mask  =ismember(matfiles,unanalyzed);
 [todo, td_order]  = sort(matfiles(todo_mask));
 ARTs           = ARTs(todo_mask);
 ARTs           = ARTs(td_order,:);
 
 % we will now create both the path list and the experiment name list
 

paths        = {};
expname      = {};
psignalfiles = {};
pathidx = 0;
expidx  = 1; 
psiidx  = 1;
for ii = 1:length(todo)
    bb = strsplit(todo{ii},'\'); % parse path 
    temppath = fullfile(file_path,bb{2}); % create path
    
    % if newpath, add it to paths and reset idices 
    if ~ismember(temppath,paths) 
        
        paths{end+1,1} = temppath;
        
        expidx  = 1;
        psiidx  = 1 ;
        pathidx = pathidx + 1; 
    end
    % add exp and psigal files 
    expname{expidx,pathidx} = strjoin(bb(3:end), '\');
    psignalfiles{psiidx,pathidx} = ARTs{ii};
    
    % increment indexes
    psiidx = psiidx + 1 ;
    expidx = expidx + 1 ;
    
end

[~,animalIDs] = cellfun(@fileparts,paths,'Uni',0);
        
    
    
   
    
    
 
     

