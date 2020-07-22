function [ActiveFolders,PsignalFiles] = CreateAnimalActiveList(AnimalID,savepath)
% this function looks to see what files are ready to be proccessed and
% extracted via the TwoPhotonPipeline's ExtractFluorescence functionality
% and looks for experiments that contain active behavior by looking at the
% naming convention of the containing folders
%
if ~exist('savepath','var')
    savepath = '//vault3/data/kelson/analyzed';
end 

if ~exist('AnimalID','var')
    AnimalPath =  uigetdir(savepath);
    AnimalID = fileparts(AnimalPath);
else 
  AnimalPath =  fullfile(savepath,AnimalID);
end

% this will currently break if you put in a nonexistant animalID
 ActiveFolders = dir( [AnimalPath, '\**\TonesNoiseActive\Fluorescence.mat']);

 if isempty(ActiveFolders) 
    % if running on data on the source psignalfiles on google drive 
    PsignalFiles = PsignalFileCheck(AnimalPath,true);    
    [ActiveFolders{1:length(PsignalFiles)}] = deal(AnimalPath);
 else
     % if running on vault data
     ActiveFolders = struct2cell(ActiveFolders);
     ActiveFolders = ActiveFolders(2,:)';
   
     ActiveFolders2 = dir( [AnimalPath, '\**\ToneNoiseActive\Fluorescence.mat']);
     ActiveFolders2 = struct2cell(ActiveFolders2);
     ActiveFolders2 = ActiveFolders2(2,:)';
     
     ActiveFolders = cat(1,ActiveFolders,ActiveFolders2);
     PsignalFiles = cellfun(@PsignalFileCheck,ActiveFolders,'UniformOutput',0);
end 


end 

function  psignal_files = PsignalFileCheck(Directory,all_matches) 

if ~exist('all_matches','var')
    all_matches = false ;
end 
PathFiles =  dir(Directory);
PathFiles = struct2cell(PathFiles);
PathFiles = PathFiles(1,:);
% find files
expression = 'ART.*mat';
if all_matches
    matches =  regexp( PathFiles,expression);
else 
   matches =  regexp( PathFiles,expression,'once');
end 

% convert logicnl index to row number
psignal_file_idx = find(~cellfun(@isempty,matches));

% extract file 
if length(psignal_file_idx) == 1 
    psignal_files = PathFiles{psignal_file_idx};
else 
    psignal_files = PathFiles(psignal_file_idx);
end 


end 


