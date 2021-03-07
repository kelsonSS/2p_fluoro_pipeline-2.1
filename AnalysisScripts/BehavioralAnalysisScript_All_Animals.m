% x = dir('C:\Users\KanoldLab\Google Drive\PsignalData\KSS\Analyzed');
% 
% x = {x.name}';

 x = cellfun(@(x) strsplit(x,'\'),TNBehavior(:,3),'UniformOutput',0);
 
 x = unique(cellfun(@(x) x(4),x));


for ii = 1:length(x)
    BehavioralAnalysisByAnimal(x{ii},'Totals');
    close all 
end 

