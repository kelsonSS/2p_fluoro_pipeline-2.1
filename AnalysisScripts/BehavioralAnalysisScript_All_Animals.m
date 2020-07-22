x = dir('C:\Users\KanoldLab\Google Drive\PsignalData\KSS\');

x = {x.name}';

for ii = 1:length(x)
    BehavioralAnalysisByAnimal(x{ii},'Totals');
    close all 
end 

