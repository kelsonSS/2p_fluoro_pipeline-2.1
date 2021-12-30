function results = ExtractAnimalBehaviorDaily(psi,daily_plot_flag)

%% plot 1: hit rate graph 
hitrate =  [psi.Performance.HitRate];
hitrate = hitrate(1:end-1);


missrate = [psi.Performance.MissRate];
missrate = missrate(1:end-1);


earlyrate = [psi.Performance.EarlyRate];
earlyrate = earlyrate(1:end-1);

earlyrate2 = [psi.Performance(1:end-1).EarlyTrial] == 3;

earlyrate2 = sum(earlyrate2) / length(earlyrate2); 
hitrate2 = hitrate(end) - earlyrate2;
dprime2 = norminv(hitrate(end)) - norminv(earlyrate2);
dprime = norminv(hitrate2(end)) - norminv(earlyrate2);


%% Plot 2 lick histogram

TrialDurSeconds = psi.PrimaryDuration + psi.PreStimSilence + psi.PostStimSilence ;

LickLatency = [psi.FirstResponse];

%% Plot 3 Performance by level 

uLevels = sort(unique(psi.Levels),'descend');
uLevelsSNR = uLevels - 50;
LevelsSNR = psi.Levels - 50;

day = 1;
for lvl = 1:length(uLevels)
    lvl_idx = psi.Levels  == uLevels(lvl);
    % find relevant trials 
    hits_lvl_temp =  psi.Hits(lvl_idx);
    early_lvl_temp = psi.Early(lvl_idx);
    num_trials(lvl) = length(hits_lvl_temp);
    percent_correct(lvl) = nansum(hits_lvl_temp) / num_trials(lvl);
    percent_early (lvl) = nansum(early_lvl_temp) / num_trials(lvl);
    dprime_level(lvl) = norminv(percent_correct(lvl)) - norminv(percent_early(lvl));
    dprime_level(lvl) = clip(dprime_level(lvl),-3,3);
end 
    percent_early_total = psi.Performance(end-1).EarlyRate;
    percent_correct_total = psi.Performance(end-1).HitRate - percent_early_total;
    percent_miss_total = psi.Performance(end-1).MissRate;
    assert ( (percent_early_total + percent_correct_total + percent_miss_total -1) < 1e-6 )
    dprime_total = norminv(percent_correct_total) - norminv(percent_early_total);
    dprime_total = clip(dprime_total,-3,3);
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
    title('Trials Correct vs. SNR')
    ylim([0 1])
    xlabel(' dB SNR')
    ylabel('Fraction Correct')


end 

results.HitRate = hitrate;
results.MissRate = missrate;
results.EarlyRate = earlyrate;
results.EarlyRate2 = earlyrate2;
results.Dprime2 = dprime2;
results.Levels = uLevels;
results.SNR = uLevelsSNR;
results.PercentCorrect = percent_correct;
results.DprimeLevel = dprime_level;
results.DprimeTotal =  dprime_total;
results.TrialsPerLevel = num_trials;
results.LickLatency = LickLatency(:,1);
results.LevelsSNR = LevelsSNR;
