% for animal 1:e
% 
% get animal info file 
% [animalName,...] = AnimalInfoParse()
% 
% Find all Active animal trials
%  [] =  FindActiveTrials()
 psi = WF_getPsignalInfo();
% psi= WF_getPsignalInfo(DataDirs)

%% plot 1: hit rate graph 
hitrate =  [psi.Performance.HitRate];
hitrate_day = hitrate(1:end-1);
figure; hold on
plot(hitrate_day)


missrate = [psi.Performance.MissRate];
missrate = missrate(1:end-1);
plot(missrate)

earlyrate = [psi.Performance.EarlyRate];
earlyrate = earlyrate(1:end-1);
plot(earlyrate)
title('hitrate, missrate, and falsealarm rate');
legend({'hitrate','missrate','falsealarmrate'})

%% Plot 2 lick histogram



%% Plot 3 Performance by level 

uLevels = unique(psi.Levels);
day = 1;
for lvl = 1:length(uLevels)
    lvl_idx = psi.Levels == uLevels(lvl);
    % find relevant trials 
    hits_lvl_temp =  psi.Hits(lvl_idx);
    num_trials(lvl) = length(hits_lvl_temp);
    percent_correct(lvl) = sum(hits_lvl_temp) / num_trials(lvl);
 
    % once adding days      
    % 
    % 
    %
    
end 
% plot
figure;bar(percent_correct)
xticklabels(uLevels); %change to SNR noise = 50db
% label x and y axis 

% title(animal info name)


