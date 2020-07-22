function CreateAnimalAgeDatabase(InPath)
% this function loops over all PsignalFolders containing an AnimalInfo.mat 
% file and creates a table containing the animals mean age during experiments 
% as well ast the age when it started experiments and the age when it
% eneded experiemnts
% 
% InPath - if blank or run on all animal data will return a table of all animal 
% info  
% 
% - if InPath is run on one animals' folder it will return the age info for
% that one animal 

% Kelson Shilling-Scrivo 2020
%
% 
%

if ~exist('InPath','var')
    InPath = 'C:\Users\KanoldLab\Google Drive\PsignalData\KSS';
end 

animal_dirs =  dir( [InPath, '\**\AnimalInfo.mat'] ) ;
animal_paths = {animal_dirs.folder};

for ii = 1:length(animal_paths);
    
     load(fullfile(animal_paths{ii},'AnimalInfo.mat'));
     
     AnimalID = AnimalInfo.Animal;
     if isempty(AnimalInfo.DOB)
         warning('AnimalName has no insufficent AnimalInfo')
         continue
     end 
     
     
     
     
    
    
end 









