function GroupData = BehavioralAnalysisByGroup()

AllSNRs = [30; 20; 10 ; 5 ; 0; -5; -10; -20]

behavior_folder = uigetdir('//vault3/data/kelson');


Animal_IDs = dir(behavior_folder);
Animal_IDs  = {Animal_IDs.name};
Animal_IDs = Animal_IDs(3:end);    



for animal_idx = 1:length(Animal_IDs)
    animal_ID = Animal_IDs{animal_idx}(1:end-4);
    load(fullfile(behavior_folder,Animal_IDs{animal_idx}));
    uLevels = resultsAll.Combined.SNR;
    
    
 if ~exist('GroupData','var')    
     GroupData = struct('Names',{{resultsAll.AnimalInfo.Animal}},...
     'Ages',{resultsAll.AnimalInfo.EndAgeMonths},...
     'NumTrainingDays',length(resultsAll.Daily),...
     'MeanAgeDays', resultsAll.AnimalInfo.MeanAgeDays,...
     'PsignalData',{resultsAll.PsignalData},...
     'AnimalInfo', {resultsAll.AnimalInfo},...
     'AllResults', PackageResults([],...
                                          resultsAll.Combined,...
                                          uLevels) , ...
     'LastWeekPerformance',PackageResults([] ,...
                                                   resultsAll.LastWeekPerformance,...
                                                   uLevels));
 else
%     
    
    
    GroupData.Names{end+1} = resultsAll.AnimalInfo.Animal;
    GroupData.Ages(end+1) = resultsAll.AnimalInfo.EndAgeMonths;
    GroupData.NumTrainingDays(end+1) = length(resultsAll.Daily);
    GroupData.MeanAgeDays(end+1) = resultsAll.AnimalInfo.MeanAgeDays;
    GroupData.AllResults = PackageResults(GroupData.AllResults,...
                                          resultsAll.Combined,...
                                          uLevels);
    if ~isempty(resultsAll.LastWeekPerformance)
        GroupData.LastWeekPerformance = PackageResults(...
                                           GroupData.LastWeekPerformance,...
                                           resultsAll.LastWeekPerformance,...
                                           resultsAll.LastWeekPerformance.SNR);
    

    
    
 end
 GroupData.LastWeekPerformance
 GroupData.AllResults
 end 
end
 GroupData.SNR = [30; 20; 10 ; 5 ; 0; -5; -10; -20];



%% Plotting 
PlotDetectionData(GroupData);


end 

function out = PackageResults(old, new,newSNRs,AllSNRs)
% this function merges structures with identical fields and all numeric data 
% and merges the two into the structure out 
if ~exist('AllSNRs','var')
AllSNRs = [30; 20; 10 ; 5 ; 0; -5; -10; -15 ;-20];
end 
% create SNR index

[~,snr_idx ,~ ]= intersect(AllSNRs,newSNRs);
                


% empty struct check

if isempty(new);out = old;return; end 


out = struct();

if ~isempty(old)
    fields = fieldnames(old);
else 
    fields = fieldnames(new);
end 

for f_idx = 1:length(fields)
    curr_f = fields{f_idx};
    %initialize new field with correct number of levels
    if strmatch(curr_f,'SNR')
        continue
    end 
    
    if iscell(new.(curr_f))
        new_field_spaced = cell(size(AllSNRs)) ;
    else
        new_field_spaced = nan(size(AllSNRs));
    end 
    % add entries to correct SNRs
    new_field_spaced(snr_idx) =  new.(curr_f); 
    if  isempty(old) || isempty(old.(curr_f))
        out.(curr_f) = new_field_spaced;
    else 
        out.(curr_f) = cat(2, old.(curr_f),new_field_spaced ) ;
    end  
        
end 

end 









