function results = ExtractAnimalBehaviorDaily(psi,daily_plot_flag)

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
LevelsSNR = psi.Levels - 50;

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
results.LevelsSNR = LevelsSNR;
