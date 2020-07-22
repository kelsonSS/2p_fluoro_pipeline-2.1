function  [resultsAll] =  BehavioralAnalysisByAnimal(animal_ID,toPlot)

% this function takes in the AnimalID which corresponds to a folder
% containing the psignalfiles of this animals behavioral performance. This
% function then will look at the daily and overall performance of the
% animal
% Inputs
% AnimalID - name of animal folder containing psignalfiles
% toPlot - optional limit number of figures to plot 
%          All - default
%          Totals - only the combined plots of animals overall 
%                  performance will be plotted 
%

%
% TODO : extend to be able to graph data generated from training box
Main_Psignal_Folder = 'C:\Users\KanoldLab\Google Drive\PsignalData\KSS';

daily_plot_flag = 1;
weekly_plot_flag = 1;
resultsAll ={};
%% input
if ~exist('toPlot','var');toPlot = 'All';end
if strcmp(toPlot,'Totals');daily_plot_flag = 0; end 

if ~exist('animal_ID','var')|| isempty(animal_ID)
    folder_path = uigetdir(Main_Psignal_Folder);
    [save_path, animal_ID] = fileparts(folder_path);
    
    [ActiveFolders, PsignalFiles] = CreateAnimalActiveList(animal_ID,save_path);
else 
 [ActiveFolders, PsignalFiles] = CreateAnimalActiveList(animal_ID,Main_Psignal_Folder);
end 

fprintf('Analyzing %s \n', animal_ID);

if ischar(PsignalFiles) || isempty(PsignalFiles)
    return
end 
%% Daily processing 
psigAll ={};
dailyResults = {};
for psig_file_idx = 1:length(ActiveFolders)  
    % get psignal info
    psig_file = fullfile(ActiveFolders(psig_file_idx),...
                         PsignalFiles(psig_file_idx));
       
     
    psigDay = WF_getPsignalInfo(char(psig_file) );
   
    % check for correct sound class
    if  ~contains(psigDay.Class,'Tone')
        continue
    end 
    % plot 
    if daily_plot_flag
        figure
        title(PsignalFiles(psig_file_idx), 'Interpreter','None')
    end 
    if isfield(psigDay,'Performance') % only look at behavior expts.
        dailyResults{end+1} = PlotAnimalBehaviorDaily(psigDay,daily_plot_flag);
     
        psigAll{end+1} = psigDay; 
    end 
end 

% run combined stats
figure 
title(sprintf( '%s Combined Stats',animal_ID))

if length(psigAll) == 0 
    return
end
%% Combined Results
combinedResults =  PlotAnimalBehaviorCombined(dailyResults,psigAll);

numTrainingDays = length(dailyResults);
if numTrainingDays >= 10  % if there were at least 2 weeks of training  
    LastFiveResults =  PlotAnimalBehaviorCombined(dailyResults(end-4:end),psigAll(end-4:end));
else 
     LastFiveResults = [];
end 
%% animal Info Parsing
% Attempt to grab AnimalInfoFile Which will load a struct called AnimalInfo
 animal_info_file = fullfile(ActiveFolders(1),'AnimalInfo.mat');
    try
       load(char(animal_info_file));
       
       load('\\vault3\data\kelson\AgeTable.mat')
          
       try
           % if animal age is in table add it 
           assert(~isnat(AgeTable{animal_ID,'DOB'}))
           AnimalInfo.DOB = AgeTable{animal_ID,'DOB'};
           
       catch
           
           
           if ~isempty(AnimalInfo.DOB) %if there is a DOB
               month_DOB = str2num(AnimalInfo.DOB(1:2));
           else
               month_DOB = -1; % bad DOB file 
           end
           if month_DOB <1 % if bad DOB
               AnimalInfo = [];
               warning('Animal info DOB incorrect for %s',animal_ID)
           elseif month_DOB < 13 % if month  presented first 
               AnimalInfo.DOB = datetime(AnimalInfo.DOB,'InputFormat', 'MM/dd/yy');
           else % day presented first 
               AnimalInfo.DOB = datetime(AnimalInfo.DOB,'InputFormat', 'dd/MM/yy');
           end
           
            % add to AgeTable 
            
            
            AgeTable{animal_ID,1} = AnimalInfo.DOB;
            
            save('//vault3/data/kelson/AgeTable.mat','AgeTable')
                 
           
           
       end  
          
           
       catch
           
           warning('AnimalInfoFile not found')
           AnimalInfo = [];
    end
       
    
           if ~isempty(AnimalInfo) 
               AnimalInfo = AnimalAgeAnalysis(AnimalInfo,psigAll{1}.ExperimentStart,...
                   psigAll{end}.ExperimentStart);
           end
    

    
    
%% packing 
resultsAll.Daily = dailyResults;
resultsAll.Combined = combinedResults;
resultsAll.PsignalData = psigAll;
resultsAll.AnimalInfo = AnimalInfo;
resultsAll.LastWeekPerformance = LastFiveResults; 
out_folder = 'Z:/TNDetectionAnalysis';

out_path = fullfile(out_folder,animal_ID);

save(out_path,'resultsAll','-v7.3')

end 

function resultsAll = PlotAnimalBehaviorCombined(results,psigAll)

%% Plot Average number of hits by Day

%hits_daily = cellfun(@(x) nanmedian(x.Hits),results,'UniformOutput',1);
%histogram(hits_daily)
%% Plot Average Latency by Day

median_lick_latency = cellfun(@(x) nanmedian(x.LickLatency),results,'UniformOutput',1);

figure 
bar(median_lick_latency)
title('Median Lick Latency')
xlabel('Day')
ylabel('Time (s)')



%% plot overall percentage of hits by level 
% find all unique Levels
ulevels_total = unique(cell2mat(cellfun(@(x) x.Levels',results,'UniformOutput',0)));
ulevels_SNR = ulevels_total - 50; 
ul_len = length(ulevels_total);
trials_lvl_total =zeros(ul_len,1);
percent_correct_total = zeros(ul_len,1);
early_rate_total = zeros(ul_len,1);
lick_latency_total = cell(ul_len,1);
for day_idx = 1:length(results)
    
 psig = results{day_idx};
    
 ulevels = psig.Levels;
 % construct percent correct and total_trials
 for psig_lvl = 1:length(ulevels)
    % find which level to add to total 
     total_lvl_idx = find(ulevels(psig_lvl) == ulevels_total);
     % extract info
     num_trials_lvl =  psig.TrialsPerLevel(psig_lvl);
     pc_lvl = psig.PercentCorrect(psig_lvl) ;
     % add to percent correct
     percent_correct_total(total_lvl_idx)= percent_correct_total(total_lvl_idx)+...
         pc_lvl * num_trials_lvl ;
     % add to total trials 
     trials_lvl_total(total_lvl_idx) = trials_lvl_total(total_lvl_idx) + num_trials_lvl;
       
     % extract early rate 
     early_rate = sum( psigAll{day_idx}.EarlyLevels{psig_lvl});
     early_rate_total(total_lvl_idx) = early_rate_total(total_lvl_idx) + ...
                                        early_rate ;
     
                                    
     % extract Lick histogram by level  
       trial_idx = psigAll{day_idx}.Levels ==  ulevels(psig_lvl);
       lick_timing = psig.LickLatency(trial_idx,1);
       lick_latency_total{total_lvl_idx} = cat(1,...
                         lick_latency_total{total_lvl_idx},lick_timing);
      
     
     % lick_histgram_levels{lvl} = cat(1 (or 2?) , lick_hisgtram{lvl} + new_licks ) 
 end

% lick_histgram = [1x 14] [1x 50] [1x 43] 
%  for ii = 1 : length(lick_histgram)
%       figure 
%          histogram(lick_histgram{ii}
%            title()
% end 
% 

end

 percent_correct_total = percent_correct_total./trials_lvl_total ;
 early_rate_total = early_rate_total./trials_lvl_total;
 

 % standard error of the mean = sqrt( (P * (1-P)) / N)
SEM_hits =  sqrt( percent_correct_total .* (1-percent_correct_total) ./ ... 
                              trials_lvl_total );                          
SEM_early = sqrt( early_rate_total .* (1-early_rate_total) ./ ... 
                              trials_lvl_total );
 
 %% plot showing fraction correct of licks by level                         
 
 % munge data into vectors x and g corresponding to lick times and lick group 
 x =[];
 g = [];
 for ii = 1:length(lick_latency_total)
     x = cat(1,x, lick_latency_total{ii});
     L = length(lick_latency_total{ii});
     g = cat(1,g, ii * ones(L,1));
 end 
 
 figure 
 boxplot(x,g)
 title('Latency by Level')
 xlabel('dB SNR')
 ylabel('Time (s)')
 xticks(1:length(ulevels_total))
 xticklabels(ulevels_SNR)
 
                          
 %% plot hits                          
 figure
 bar(percent_correct_total)
 hold on
 xticks(1:length(ulevels_total))
 xticklabels(ulevels_SNR);
 errorbar(percent_correct_total,SEM_hits ,'.')
 title('Fraction Correct vs. dB SNR')
 ylabel('Fraction Correct')
 xlabel ('dB SNR')
 
 %% plot for total number of trials per level 
 figure
 bar (trials_lvl_total)
 xticklabels(ulevels_SNR);
 title('Number of Trials per dB SNR')
 ylabel('Number of Trials')
 xlabel('dB SNR')
 
 %% plot of error rate by trial
 figure
 bar(early_rate_total)
 hold on
 errorbar(early_rate_total,SEM_early ,'.')
 title('Error Rate by dB SNR')
 xticklabels(ulevels_SNR);
 xlabel('dB SNR')
 ylabel('Fraction')
 
 
 
 
 % plot D'prime by level 
%  norminv(HitRate) - norminv(EarlyRate)
  
  dprime= norminv(percent_correct_total) - norminv(early_rate_total);

 figure
 bar(dprime)
 title('Sensitivity')
 xticklabels(ulevels_SNR);
 ylabel('D prime')
 xlabel('dB SNR')
% make function to get earlyRate by Level 


resultsAll.SNR = ulevels_SNR;
resultsAll.DPrime = dprime;
resultsAll.HitRateMean = percent_correct_total;
resultsAll.HitsRateSEM = SEM_hits;
resultsAll.EarlyRateMean = early_rate_total;
resultsAll.EarlyRateSEM = SEM_early;
resultsAll.NumTrials =trials_lvl_total;
resultsAll.LickLatency = lick_latency_total;



end 

    function results = PlotAnimalBehaviorDaily(psi,daily_plot_flag)

%% plot 1: hit rate graph 
hitrate =  [psi.Performance.HitRate];
hitrate = hitrate(1:end-1);

missrate = [psi.Performance.MissRate];
missrate = missrate(1:end-1);


earlyrate = [psi.Performance.EarlyRate];
earlyrate = earlyrate(1:end-1);


%% Plot 2 lick histogram

TrialDurSeconds = psi.PrimaryDuration + psi.PreStimSilence + psi.PostStimSilence ;

LickLatency = [psi.FirstResponse];

%% Plot 3 Performance by level 

uLevels = unique(psi.Levels);
uLevelsSNR = uLevels - 50;

day = 1;
for lvl = 1:length(uLevels)
    lvl_idx = psi.Levels  == uLevels(lvl);
    % find relevant trials 
    hits_lvl_temp =  psi.Hits(lvl_idx);
    num_trials(lvl) = length(hits_lvl_temp);
    percent_correct(lvl) = nansum(hits_lvl_temp) / num_trials(lvl);
    
end 
% plot
if daily_plot_flag
     % fig 1     
    h1= figure; hold on
    plot(hitrate)
    plot(missrate)
    plot(earlyrate)
    title('hitrate, missrate, and falsealarm rate');
    legend({'hitrate','missrate','falsealarmrate'})
    hold off;
    % fig 2
    h2 = figure; histogram(LickLatency,'BinWidth',.1)
    xlim([0 TrialDurSeconds])
    title('Lick Rate over Trial Duration')
    xlabel('Time in Seconds')
    ylabel('Licks')

    % fig 3 
    h3 = figure; hold  on; bar(percent_correct)
    xticks(1:length(uLevels))
    xticklabels( num2str(uLevelsSNR) ); %change to SNR noise = 50db
    title('Trial Correctness vs. SNR')
    ylim([0 1])
    xlabel(' dB SNR')
    ylabel('Fraction Correct')


end 

results.HitRate = hitrate;
results.MissRate = missrate;
results.EarlyRate = earlyrate;
results.Levels = uLevels;
results.LevelsInSNR = uLevelsSNR;
results.PercentCorrect = percent_correct;
results.TrialsPerLevel = num_trials;
results.LickLatency = LickLatency(:,1);

end 
