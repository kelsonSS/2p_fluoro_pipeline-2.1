function Info = getAnimalInfo(AnimalList,BasePath)
% this function gives extracts animal info from animals used during experiments 
% This function creates a table containing the animals age during experiments 
% as well as sex and DOB
% 
% AnimalIDs - the folder names for each animal where all the psignalfiles
% are stored
% 
% BasePath - Path to the users PsignalData Folder

% Kelson Shilling-Scrivo 2020
%
% 
%

AnimalID = {};
DOB = {};
Sex = {};
Experiment_Date ={};
Age_Months = {};
Age_Days = [];

if ~exist('BasePath','var')
    BasePath= 'C:\Users\kelsonss\Google Drive\psignalData\KSS';
end 



for ii = 1:length(AnimalList)
    
    AnimalName =  AnimalIDFromPath(AnimalList{ii});
    
    if iscell(AnimalName)
        AnimalName = AnimalName{1};
    end 
        
    try  
    % loading 
     load(fullfile(BasePath,AnimalName,'AnimalInfo.mat'));
     PsignalFile = PsignalFileCheck(AnimalList{ii});
     handles =WF_getPsignalInfo(fullfile(AnimalList{ii},PsignalFile));    
        
     % animal ID
     AnimalID{ii} = AnimalInfo.Animal;
     
     % Age
         Expt_date = handles.ExperimentDate;
     
     try
     DOB{ii} = datetime(AnimalInfo.DOB , 'InputFormat' ,'dd/MM/yy');
     
     % we never use animals below 40  days old
     assert(day(caldiff([DOB{ii},Expt_date],'days')) > 60)

     
     catch 
          DOB{ii} = datetime(AnimalInfo.DOB , 'InputFormat' ,'MM/dd/yy');
     end 
         
     Expt_date = handles.ExperimentDate;
     
     Experiment_Date{ii} = Expt_date;
     Age_Months{ii} = between(DOB{ii},Expt_date,{'months','days'});
     Age_Days(ii) = caldays(caldiff([DOB{ii},Expt_date],'days'));
     

     
     % Sex
     Sex{ii} = upper(AnimalInfo.sex);
     
    catch 
         warning(sprintf('%s has no insufficent AnimalInfo',AnimalList{ii}))
         continue
     end 
     
end 
     
    Info = table(AnimalID',Sex',DOB',Experiment_Date',Age_Months',Age_Days','VariableNames',...
                 {'AnimalID','Sex','DOB','ExptDate','Months_old','Days_old'}); 
    
    
end 









