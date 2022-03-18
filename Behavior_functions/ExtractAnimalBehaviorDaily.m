function results = ExtractAnimalBehaviorDaily(psi,daily_plot_flag)

%% plot 1: hit rate graph 


missrate_over_time = [psi.Performance.MissRate];
missrate_over_time = missrate_over_time(1:end-1);

earlyrate_over_time = [psi.Performance.EarlyRate];
earlyrate_over_time = earlyrate_over_time(1:end-1);

hitrate_over_time = [psi.Performance.HitRate];
hitrate_over_time = hitrate_over_time(1:end-1);



%% Plot 2 lick histogram

TrialDurSeconds = psi.PrimaryDuration + psi.PreStimSilence + psi.PostStimSilence ;

LickLatency = [psi.FirstResponse];

%% Plot 3 Performance by level 
num_trials_total = psi.Performance(end).Trials;
uLevels = sort(unique(psi.Levels),'descend');
uLevelsSNR = uLevels - 50;
LevelsSNR = psi.Levels - 50;

day = 1;
for lvl = 1:length(uLevels)
    lvl_idx = psi.Levels  == uLevels(lvl);
    % find relevant trials 
    hits_lvl_temp =  psi.Hits(lvl_idx);
    early_lvl_temp = psi.Early(lvl_idx);
    num_trials(lvl) = sum(lvl_idx);
    percent_correct(lvl) = nansum(hits_lvl_temp) / num_trials(lvl);
    percent_early (lvl) = nansum(early_lvl_temp) / num_trials(lvl);
    dprime_level(lvl) = getCorrectedDPrime(...
                               percent_correct(lvl),...
                                percent_early(lvl),...
                                num_trials(lvl));
end 
    percent_early_total = sum(psi.Early) / num_trials_total;
    percent_correct_total =sum(psi.Hits) / num_trials_total;
    percent_miss_total = sum(psi.Miss) / num_trials_total;
    assert ( (percent_early_total + percent_correct_total + percent_miss_total -1) < 1e-6 )
    dprime_total = getCorrectedDPrime(percent_correct_total,percent_early_total,num_trials_total);
    
% plot
if daily_plot_flag
     % fig 1     
    h1= figure; hold on
    plot(hitrate_over_time)
    plot(missrate_over_time)
    plot(earlyrate_over_time)
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

results.HitRate = hitrate_over_time;
results.MissRate = missrate_over_time;
results.EarlyRate = earlyrate_over_time;
results.Levels = uLevels;
results.SNR = uLevelsSNR;
results.PercentCorrect = percent_correct;
results.PercentEarly = percent_early;
results.DprimeLevel = dprime_level;
results.DprimeTotal = dprime_total;
results.TrialsPerLevel = num_trials;
results.LickLatency = LickLatency(:,1);
results.LevelsSNR = LevelsSNR;




function dPrime = getCorrectedDPrime(hitrate,earlyrate,num_trials)
% this function extracts the hit, miss or early rate using the 
% log-linear transform for extremes ( no hits or all hits)
% 
% essentially, we add half a trial to each condition such that:
% lowest possible performance  = .5/(num_trials +1)
% best performance  = 1 - ( .5/(num_trials +1))
%
% Hautus, Micheal (1995) Corrections for Extreme Proportions... 
% Behavioral Research Methods 

hitrate_C = ((hitrate * num_trials) + .5)/ (num_trials+1);
earlyrate_C = ((earlyrate * num_trials) + .5)/ (num_trials+1);

dPrime = norminv(hitrate_C) - norminv(earlyrate_C);
dPrime = min(dPrime,2.9);
assert(~ isnan(dPrime))










