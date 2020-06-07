function NeedsPackaging(staging_folder_path)
% this function looks for files in the current staging area and returns a
% list of folders that are currently in that staging area. In this model
% Files should all be moved to the ready for processing area. the current
% model is below
%
% Staging area = Transfer folder 
%              |
%              V
% Ready to process = Files-to-unload
%              |
%              V
%      analyzed = Analyzed
%

if ~exist('staging_folder_path','var')
    staging_folder_path = '\\vault3\data\Transferfolder\forKelson' 
    x = dir(' 
