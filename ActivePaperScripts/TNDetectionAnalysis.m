% BehaviorAnalysis Script

[BehaviorListAll,UnPaired] = FindPairedBehaviorExperiments;

TN_index = contains(BehaviorListAll(:,1),'Noise');
TNBehaviorList = BehaviorListAll(TN_index,:); 
TNBehavior = Fluoro_to_Table_Behavior(TNBehaviorList);

% remove any rows that didn't properly get extracted 
is_paired_idx = ~cellfun(@isempty,TNBehavior(:,1)) &...
                ~cellfun(@isempty,TNBehavior(:,2)) ;
            
TNBehavior = TNBehavior(is_paired_idx,:);

%  ensure experiment dates are after 2019

%TNBehavior = TNBehavior( cellfun(@(x) x.Year >= 2020,...
%                                  TNAnimalInfo.ExptDate),:); 

% Behavior has some files with different amount of PreStim Silence
% trim them all to all be consistent with 1 second prestim silence 
TNBehavior(:,1) = TNBehavior_Standardize_Freqs(TNBehavior(:,1));
TNBehavior(:,2) =  TNBehavior_Fix_Timings(TNBehavior(:,2));
TNBehaviorInfo = getAnimalInfo(TNBehavior(:,4));

%% Create full BehaviorTable
TNFullBehaviorTable = CreateFullAnimalExperimentSheet(TNBehavior)

% move licks to behavior
TNBehavior(:,2) = AddLicks_to_Behavior(TNBehavior(:,2),TNFullBehaviorTable.LickResponseTimes);

%% Create Various subsets
dPrime_gt1_idx = cellfun(@(x) any(x >= 1), TNFullBehaviorTable.dPrimeLevel);
old_idx = contains(lower(TNBehavior(:,3)),'ia');
male_idx = strcmp(TNBehaviorInfo.Sex,'M');
female_idx = strcmp(TNBehaviorInfo.Sex,'F');


% subset into different groups
TNBehaviorYoung = TNBehavior(~old_idx,:);
TNBehaviorYoung_Male = TNBehavior((~old_idx) & male_idx ,:);
TNBehaviorYoung_Female = TNBehavior((~old_idx) & (~ male_idx) ,:);

TNBehaviorOld = TNBehavior(old_idx,:);
TNBehaviorOld_Male = TNBehavior(old_idx & male_idx ,:);
TNBehaviorOld_Female = TNBehavior(old_idx & (~ male_idx) ,:);





TNAgingAnimalInfo = getAnimalInfo(TNBehaviorOld(:,4));
TNYoungAnimalInfo = getAnimalInfo(TNBehaviorYoung(:,4));


%% create output struct
BehaviorResults = struct();

BehaviorResults.Young.CellInfo = CreateBehaviorCellList(TNBehaviorYoung);
BehaviorResults.Old.CellInfo = CreateBehaviorCellList(TNBehaviorOld);

AnimalIDS = cellfun(@(x) strsplit(x,filesep),TNBehavior(:,3),'UniformOutput',0);
 AnimalIDS = unique( cellfun(@(x) x{4},AnimalIDS,'UniformOutput',0));
 BehaviorResults.AnimalID = AnimalIDS;
 clear AnimalIDS
 
%% clustering 
% BehaviorResults.Young.Clust.Passive = Cluster_DF(TNBehaviorYoung(:,1),'K-means');
% BehaviorResults.Young.Clust.Active = Cluster_DF(TNBehaviorYoung(:,2),'K-means',9,...
% 'normalized',BehaviorResults.Young.Clust.Passive.Centroids);
% 
% BehaviorResults.Old.Clust.Passive = Cluster_DF(TNBehaviorOld(:,1),'K-means');
% BehaviorResults.Old.Clust.Active = Cluster_DF(TNBehaviorOld(:,2),'K-means',9,...
% 'normalized',BehaviorResults.Old.Clust.Passive.Centroids);
% 
% BehaviorResults.Combined.Clust.Passive = Cluster_DF(TNBehavior(:,1),'K-means',10);
% BehaviorResults.Combined.Clust.Active = Cluster_DF(TNBehavior(:,2),'K-means',6,...
% 'normalized',BehaviorResults.Combined.Clust.Passive.Centroids);


%% plot Cluster transitions 
% Young_passive_clusters = BehaviorResults.Young.Clust.Passive.Clusters(...
%      BehaviorResults.Young.Clust.Passive.Neuron);
%  Young_active_clusters = BehaviorResults.Young.Clust.Active.Clusters(...
%      BehaviorResults.Young.Clust.Active.Neuron);
%  
%  Old_passive_clusters = BehaviorResults.Old.Clust.Passive.Clusters(...
%      BehaviorResults.Old.Clust.Passive.Neuron);
%  Old_active_clusters = BehaviorResults.Old.Clust.Active.Clusters(...
%      BehaviorResults.Old.Clust.Active.Neuron);
% 
%  % plot Young
%  trans_mat_young = histcounts2(Young_passive_clusters,Young_active_clusters(1:2892));
% trans_mat_young = trans_mat_young ./ sum(trans_mat_young(:)) ;
%  figure; imagesc(trans_mat_young)
 % Plot Old 
 
 %%
 
  behavior_path ='Z:\Kelson\TNDetectionAnalysis';
% load in behavior 

for ii = 1:length(BehaviorResults.AnimalID)
    curr_id = BehaviorResults.AnimalID{ii};
        
    animal_sex = unique(TNBehaviorInfo{strcmp(TNBehaviorInfo.AnimalID,lower(curr_id)),2})
    


    % if animal is not in path, add it to path
    if  exist(fullfile(behavior_path,[curr_id '.mat']),'file')
           BehavioralAnalysisByAnimal(curr_id,'Totals');
          % pause
    end 
        
     

    if contains(curr_id,'IA','IgnoreCase',true)
        
        
        
        BehaviorResults.Old.Behavior.All.( sprintf('animal_%s',curr_id) ) = ...
        load(fullfile(behavior_path,curr_id));
        
        BehaviorResults.Old.Behavior.(animal_sex{:}).( sprintf('animal_%s',curr_id) ) = ...
        load(fullfile(behavior_path,curr_id));


    
              
    else 
         BehaviorResults.Young.Behavior.All.( sprintf('animal_%s',curr_id) ) = ...
        load(fullfile(behavior_path,curr_id));

          BehaviorResults.Young.Behavior.(animal_sex{:}).( sprintf('animal_%s',curr_id) ) = ...
        load(fullfile(behavior_path,curr_id));
    end 
        
    
        
end 

TN_passive =Munge_DF(TNBehavior(:,1));

% analyze behavior
TNBehaviorResultsYoung = MungeBehaviorGroupData(BehaviorResults.Young.Behavior.All);
TNBehaviorResultsOld = MungeBehaviorGroupData(BehaviorResults.Old.Behavior.All);

BehaviorResults.Young.Behavior.All.group = TNBehaviorResultsYoung;
BehaviorResults.Old.Behavior.All.group = TNBehaviorResultsOld;

TNBehaviorResultsYoungMale = MungeBehaviorGroupData(BehaviorResults.Young.Behavior.M);
TNBehaviorResultsOldMale = MungeBehaviorGroupData(BehaviorResults.Old.Behavior.M);

TNBehaviorResultsYoungFemale = MungeBehaviorGroupData(BehaviorResults.Young.Behavior.F);
TNBehaviorResultsOldFemale = MungeBehaviorGroupData(BehaviorResults.Old.Behavior.F);


BehaviorResults.Young.Behavior.All.group = TNBehaviorResultsYoung;
BehaviorResults.Old.Behavior.All.group = TNBehaviorResultsOld;



%Plot group Behavior 
PlotGroupedData(BehaviorResults.Old.Behavior.group)
PlotGroupedData(BehaviorResults.Young.Behavior.group)


% analyse Passive Data
TN_young_passive = Fluoro_to_Table(TNBehaviorYoung(:,3),50000);
TN_young_passive.AnimalInfo = getAnimalInfo(TN_young_passive.DataDirs);
TN_old_passive =  Fluoro_to_Table(TNBehaviorOld(:,3),50000);
TN_old_passive.AnimalInfo = getAnimalInfo(TN_old_passive.DataDirs);

TNYoungBehaviorTable = CreateFullAnimalExperimentSheet(TNBehaviorYoung);
TNOldBehaviorTable = CreateFullAnimalExperimentSheet(TNBehaviorOld);

 % generate datasets seperated by detection prob.
 young_good_p_idx = TNYoungBehaviorTable.DprimeTotal >= 1;
 old_good_p_idx = TNOldBehaviorTable.DprimeTotal >= 1;
 old_good_p_idx(3) = 0;
 
 TN_young_good_p = Fluoro_to_Table( TN_young_passive.DataDirs(young_good_p_idx),50000);
 TN_young_bad_p = Fluoro_to_Table( TN_young_passive.DataDirs(~young_good_p_idx),50000);
 TN_old_good_p = Fluoro_to_Table(TN_old_passive.DataDirs(old_good_p_idx),50000);
 TN_old_bad_p = Fluoro_to_Table( TN_old_passive.DataDirs(~old_good_p_idx),50000);
  
 % generate datasets seperated by sex
  old_F_idx = contains(TN_old_passive.AnimalInfo{:,2},'F'); 
  young_F_idx= contains(TN_young_passive.AnimalInfo{:,2},'F');
 
  
  TN_old_passive_F = Fluoro_to_Table(TN_old_passive.DataDirs(old_F_idx),50000);
   TN_old_passive_M = Fluoro_to_Table(TN_old_passive.DataDirs(~old_F_idx),50000);
% 
%   
    TN_young_passive_F = Fluoro_to_Table(TN_young_passive.DataDirs(young_F_idx),50000);
   TN_young_passive_M = Fluoro_to_Table(TN_young_passive.DataDirs(~young_F_idx),50000);

% quick patch for youngF FLO

if size(TN_young_passive_F.FreqLevelOrder{:,1}== 80)
    idx = mod((1:80)-1,10) >= 5;
    TN_young_passive_F.FreqLevelOrder = TN_young_passive_F.FreqLevelOrder(idx,:);
    clear idx
end


if size(TN_young_passive.FreqLevelOrder{:,1}== 80)
    idx = mod((1:80)-1,10) >= 5;
    TN_young_passive.FreqLevelOrder = TN_young_passive.FreqLevelOrder(idx,:)
    clear idx
end 

%% Detection Analysis/ d Prime Analysis
detect = DetectionAnalysis(TNBehavior(:,2));
% get plots by age
DetectionAnalysis(TNBehavior(:,2),old_idx & dPrime_gt1_idx)
title('Aging Animals') 

DetectionAnalysis(TNBehavior(:,2),~old_idx)
title('Young Animals')

figure;plot(fitlm(TNFullBehaviorTable(~old_idx &dPrime_gt1_idx  ,:),'mean_single_neuron_prediction~DprimeTotal'))
xlim([0 3])
ylim([0 .7])
saveas(gcf,'MeanSingleNeuronDetectionModel_Dprime_Young.pdf')


figure;plot(fitlm(TNFullBehaviorTable(old_idx ,:),'mean_single_neuron_prediction~DprimeTotal'))
xlim([0 3])
ylim([0 .7])
saveas(gcf,'MeanSingleNeuronDetectionModel_Dprime_Aging.pdf')

dPrimeAgeModel = fitlm(TNFullBehaviorTable,'DprimeTotal~Days_old')
figure;plot(dPrimeAgeModel)
saveas(gcf,'dPrimeByAge.pdf')
figure;plot(fitlm(TNFullBehaviorTable(old_idx,:),'mean_single_neuron_prediction~DprimeTotal'))
xlim([-3 3])
ylim([0 .7])
saveas(gcf,'MeanSingleNeuronDetectionModel_Dprime_Aging.pdf')

%% Basic hit/ miss analysis.

hit_rates = {TNFullBehaviorTable{~old_idx & dPrime_gt1_idx,'HitRate'},TNFullBehaviorTable{old_idx & dPrime_gt1_idx,'HitRate'}}
miss_rates = {TNFullBehaviorTable{~old_idx & dPrime_gt1_idx,'MissRate'},TNFullBehaviorTable{old_idx & dPrime_gt1_idx,'MissRate'}}
early_rates = {TNFullBehaviorTable{~old_idx & dPrime_gt1_idx,'EarlyRate'},TNFullBehaviorTable{old_idx & dPrime_gt1_idx,'EarlyRate'}}


PlotScatterLevels({'Young-Good','Old-Good'},'All','HitRates',hit_rates{1}',hit_rates{2}')
PlotScatterLevels({'Young-Good','Old-Good'},'All','MissRates',miss_rates{1}',miss_rates{2}')
PlotScatterLevels({'Young-Good','Old-Good'},'All','EarlyRates',early_rates{1}',early_rates{2}')

%% create d' Table
DprimeTable = createDprimeTable(TNFullBehaviorTable)
DprimeTable = array2table(DprimeTable,'VariableNames',{'SNR_minus_10','SNR_minus_5','SNR_0','SNR_5','SNR_10','SNR_20'});
%% d' Analysis 
Compare2Anova(ReplaceInf(TNBehaviorResultsYoung.DPrime),...
              ReplaceInf(TNBehaviorResultsOld.DPrime),{'Young','Aging'},...
              {'20','10','0'},'DPrime')

Compare2Anova(clip(TNBehaviorResultsYoungMale.DPrime,-3,3),...
              clip(TNBehaviorResultsOldMale.DPrime,-3,3),{'Young-Male','Aging-Male'},...
              {'20','10','0'},'DPrime')
          
Compare2Anova(clip(TNBehaviorResultsYoungFemale.DPrime,-3,3),...
              clip(TNBehaviorResultsOldFemale.DPrime,-3,3),{'Young-Female','Aging-Female'},...
              {'20','10','0'},'DPrime')          
          

          



%% compare licking responses
PlotFirstResponse(TNBehavior((~old_idx) & (dPrime_gt1_idx) ,2),'Young_Good_Responses')
PlotFirstResponse(TNBehavior((~old_idx) & (~dPrime_gt1_idx) ,2),'Young_Bad_Responses')
PlotFirstResponse(TNBehavior((old_idx) & (dPrime_gt1_idx) ,2),'Old_Good_Responses')
PlotFirstResponse(TNBehavior((old_idx) & (~dPrime_gt1_idx) ,2),'Old_Bad_Responses')

%% Active Passive Gain 
Gain = struct();
Gain.YoungGood = getActivePassiveGain(TNBehavior((~old_idx) & (dPrime_gt1_idx),:),'YoungGood')
Gain.YoungBad =  getActivePassiveGain(TNBehavior((~old_idx) & (~dPrime_gt1_idx),:),'YoungBad')
Gain.OldGood =  getActivePassiveGain(TNBehavior((old_idx) & (dPrime_gt1_idx),:),'OldGood')
Gain.OldBad =  getActivePassiveGain(TNBehavior((old_idx) & (~dPrime_gt1_idx),:),'OldBad')

Compare2AnovaBehavior({Gain.YoungGood.PosCell_Gain ,Gain.YoungBad.PosCell_Gain},...
                      {Gain.OldGood.PosCell_Gain ,Gain.OldBad.PosCell_Gain},...
                      {'Young','Old'},{'Good','Bad'},'AttentionalGain_Pos')

Compare2AnovaBehavior({Gain.YoungGood.NegCell_Gain ,Gain.YoungBad.NegCell_Gain},...
                      {Gain.OldGood.NegCell_Gain ,Gain.OldBad.NegCell_Gain},...
                      {'Young','Old'},{'Good','Bad'},'AttentionalGain_Neg')                  
                  
                  


%% TemporalAnalysis

Temporal = []

Temporal.YoungGood = TemporalAnalysisByAnimal(TN_young_good_p,'all','young-good-P');
Temporal.YoungBad = TemporalAnalysisByAnimal(TN_young_bad_p,'all','young-bad-P');

Temporal.OldGood = TemporalAnalysisByAnimal(TN_old_good_p,'all','Old-good-P');
Temporal.OldBad = TemporalAnalysisByAnimal(TN_old_bad_p,'all','Old-bad-P');

% Balance Index
%% HitMiss Analysis 
getFluoroHitsMiss(TNBehavior( (dPrime_gt1_idx)& (~old_idx),2),'Young_Good')
getFluoroHitsMiss(TNBehavior( (~dPrime_gt1_idx)& (~old_idx),2),'Young_Bad')
getFluoroHitsMiss(TNBehavior( (~dPrime_gt1_idx)& (old_idx),2),'Old_Bad')
getFluoroHitsMiss(TNBehavior( (dPrime_gt1_idx)& (old_idx),2),'Old_Good')
%Intensity
 %% Intensity
   Intensity = struct()
  % Age Intensity

  % performance 
  
  FluoroIntensityAnalysis(TN_young_good_p,TN_young_bad_p,...
      'max',{'Young_Good','Young_Bad'})
  
   FluoroIntensityAnalysis(TN_old_good_p,TN_old_bad_p,...
      'max',{'Old_Good','Old_Bad'})
  
% Age-Rev
Intensity.Stats.YoungAging.Mean =  FluoroIntensityAnalysis(TN_young_passive,TN_old_passive,...
      'mean',{'Young','Aging'})
  
  
  Intensity.Stats.YoungAging.Max =  FluoroIntensityAnalysis(TN_young_passive,TN_old_passive,...
      'max',{'Young','Aging'})
   



  
  % Sex 
  
    Intensity.Stats.YoungMF.Mean = FluoroIntensityAnalysis(TN_young_passive_M,TN_young_passive_F,...
      'mean',{'Young-M','Young-F'})
  
   Intensity.Stats.YoungMF.Max =  FluoroIntensityAnalysis(TN_young_passive_M,TN_young_passive_F,...
      'max',{'Young-M','Young-F'})
  
  % old M vs young M
  
  
   Intensity.Stats.OldYoungM.Mean = FluoroIntensityAnalysis(TN_young_passive_M,TN_old_passive_M,...
      'mean',{'Young-M','Old-M'})
  
   Intensity.Stats.OldYoungM.Max =  FluoroIntensityAnalysis(TN_young_passive_M,TN_old_passive_M,...
      'max',{'Young-M','Old-M'})
  

    % old F vs young F
  
  
   Intensity.Stats.OldYoungF.Mean = FluoroIntensityAnalysis(TN_young_passive_F,TN_old_passive_F,...
      'mean',{'Young-F','Old-F'})
  
   Intensity.Stats.OldYoungF.Max =  FluoroIntensityAnalysis(TN_young_passive_F,TN_old_passive_F,...
      'max',{'Young-F','Old-F'})
  
  
  
save 'IntensityStats.mat' 'Intensity' 

 %% Timing
    %% Timing Noise 
   
    
    
  TemporalAnalysisByAnimal(TN_young_good_p,'all','Young-good-P') 
  TemporalAnalysisByAnimal(TN_young_bad_p,'all','Young-bad-P')
  TemporalAnalysisByAnimal(TN_old_good_p,'all','Old-good-P')
  TemporalAnalysisByAnimal(TN_old_bad_p,'all','Old-bad-P')
  
    
Timings = struct();
Timings.Passive.YoungF = TemporalAnalysisByAnimal(TN_young_passive_F,'all','Timing_young_F_noise');
    
Timings.Passive.YoungM = TemporalAnalysisByAnimal(TN_young_passive_M,'all','Timing_young_M_noise');

Timings.Passive.AgingF = TemporalAnalysisByAnimal(TN_old_passive_F,'all','Timing_Old_F_noise');

Timings.Passive.AgingM =TemporalAnalysisByAnimal(TN_old_passive_M,'all','Timing_Old_M_noise');


Timings.Passive.Aging = TemporalAnalysisByAnimal(TN_old_passive,'all','Timing_Old_All_noise');

Timings.Passive.Young =TemporalAnalysisByAnimal(TN_young_passive,'all','Timing_Young_All_noise');


PlotScatterLevels({'Young','Aging','YoungM','AgingM','YoungF','AgingF'},...
                   'Noise','Tone_ONOFF Bias_index_Noise',...
                Timings.Noise.Young.ToneBiasByAnimal',...
                Timings.Noise.Aging.ToneBiasByAnimal',...
                Timings.Noise.YoungM.ToneBiasByAnimal',...
                Timings.Noise.AgingM.ToneBiasByAnimal',...
                Timings.Noise.YoungF.ToneBiasByAnimal',...
                Timings.Noise.AgingF.ToneBiasByAnimal')


PlotScatterLevels({'Young','Aging','YoungM','AgingM','YoungF','AgingF'},...
                   'Noise','Noise_Tone_Onset_Bias_index',...
                Timings.Noise.Young.NoiseBiasByAnimal',...
                Timings.Noise.Aging.NoiseBiasByAnimal',...
                Timings.Noise.YoungM.NoiseBiasByAnimal',...
                Timings.Noise.AgingM.NoiseBiasByAnimal',...
                Timings.Noise.YoungF.NoiseBiasByAnimal',...
                Timings.Noise.AgingF.NoiseBiasByAnimal')


% Compare2Anova(young_timing_noise_F',young_timing_noise_M',{'Female','Male'},[1:size(young_timing_noise_M,2)],...
%     'TemporalAnalysis-MF')
% 
% 
%       % Aging
% 
% 
% Compare2Anova(old_timing_noise',young_timing_noise',{'Aging','Young'},[1:size(young_timing_noise_M,2)],...
%     'TemporalAnalysis-YoungAging')
% 
%       % Old_M vs Young_M
%       Compareare2Anova(young_timing_noise_M',old_timing_noise_M',{'Males-Young','Aging'},[1:size(old_timing_noise_M,2])
%       
% % Old_F vs Young_F
%       CompareOldYoung(young_timing_noise_F',old_timing_noise_F',[1:size(old_timing_noise_M,2)])
%       
      


save('TimingStats','Timings')


%% Bandwidth Comparision
BandwithAll = struct()

% performance 
BandwidthAnalysis(TN_young_good_p, 'interp',0,'Pos',.5,'YoungGood-Pos')
BandwidthAnalysis(TN_young_bad_p, 'interp',0,'Pos',.5,'YoungBad-Pos')
BandwidthAnalysis(TN_old_good_p, 'interp',0,'Pos',.5,'OldGood-Pos')
BandwidthAnalysis(TN_old_bad_p, 'interp',0,'Pos',.5,'OldBad-Pos')

BandwidthAnalysis(TN_young_good_p, 'interp',0,'Neg',.5,'YoungGood-Neg')
BandwidthAnalysis(TN_young_bad_p, 'interp',0,'Neg',.5,'YoungBad-Neg')
BandwidthAnalysis(TN_old_good_p, 'interp',0,'Neg',.5,'OldGood-Neg')
BandwidthAnalysis(TN_old_bad_p, 'interp',0,'Neg',.5,'OldBad-Neg')


% Using BRFS 
BandwidthAll.Young.All.Pos = BandwidthAnalysis(TN_young_passive, 'BRFS',0,'Pos',.5,'YoungAll-Pos')
BandwidthAll.Young.All.Neg = BandwidthAnalysis(TN_young_passive, 'BRFS',0,'Neg',.5,'YoungAll-Neg')

BandwidthAll.Old.All.Pos = BandwidthAnalysis(TN_old_passive, 'BRFS',0,'Pos',.5,'OldAll-Pos')
BandwidthAll.Old.All.Neg =BandwidthAnalysis(TN_old_passive, 'BRFS',0,'Neg',.5,'OldAll-Neg')

BandwidthAll.Old.All.Pos = BandwidthAnalysis(TN_old_passive, 'BRFS',0,'Pos',.5,'OldAll-Pos')
BandwidthAll.Old.All.Neg =BandwidthAnalysis(TN_old_passive, 'BRFS',0,'Neg',.5,'OldAll-Neg')



BandwidthAll.Young.F.Pos = BandwidthAnalysis(TN_young_passive_F, 'BRFS',0,'Pos',.5,'YoungF-Pos')
BandwidthAll.Young.F.Neg =BandwidthAnalysis(TN_young_passive_F, 'BRFS',0,'Neg',.5,'YoungF-Neg')

BandwidthAll.Young.M.Pos = BandwidthAnalysis(TN_young_passive_M, 'BRFS',0,'Pos',.5,'YoungM-Pos')
BandwidthAll.Young.M.Neg =BandwidthAnalysis(TN_young_passive_M, 'BRFS',0,'Neg',.5,'YoungM-Neg')


BandwidthAll.Old.F.Pos = BandwidthAnalysis(TN_old_passive_F, 'BRFS',0,'Pos',.5,'OldF-Pos')
BandwidthAll.Old.F.Neg =BandwidthAnalysis(TN_old_passive_F, 'BRFS',0,'Neg',.5,'OldF-Neg')

BandwidthAll.Old.M.Pos = BandwidthAnalysis(TN_old_passive_M, 'BRFS',0,'Pos',.5,'OldM-Pos')
BandwidthAll.Old.M.Neg =BandwidthAnalysis(TN_old_passive_M, 'BRFS',0,'Neg',.5,'OldM-Neg')


% using traditional_BW
BandwidthAll.Young.All.BW_Pos = BandwidthAnalysis(TN_young_passive, 'interp',0,'Pos',.5,'YoungAll-BW-Pos')
BandwidthAll.Young.All.BW_Neg =BandwidthAnalysis(TN_young_passive, 'interp',0,'Neg',.5,'YoungAll-BW-Neg')

BandwidthAll.Old.All.BW_Pos = BandwidthAnalysis(TN_old_passive, 'interp',0,'Pos',.5,'OldAll-BW-Pos')
BandwidthAll.Old.All.BW_Neg =BandwidthAnalysis(TN_old_passive, 'interp',0,'Neg',.5,'OldAll-BW-Neg')

BandwidthAll.Old.All.BW_Pos = BandwidthAnalysis(TN_old_passive, 'interp',0,'Pos',.5,'OldAll-BW-Pos')
BandwidthAll.Old.All.BW_Neg =BandwidthAnalysis(TN_old_passive, 'interp',0,'Neg',.5,'OldAll-BW-Neg')



BandwidthAll.Young.F.BW_Pos = BandwidthAnalysis(TN_young_passive_F, 'interp',0,'Pos',.5,'YoungF-BW_Pos')
BandwidthAll.Young.F.BW_Neg =BandwidthAnalysis(TN_young_passive_F, 'interp',0,'Neg',.5,'YoungF-BW_Neg')

BandwidthAll.Young.M.BW_Pos = BandwidthAnalysis(TN_young_passive_M, 'interp',0,'Pos',.5,'YoungM-BW_Pos')
BandwidthAll.Young.M.BW_Neg =BandwidthAnalysis(TN_young_passive_M, 'interp',0,'Neg',.5,'YoungM-BW_Neg')

BandwidthAll.Old.M.BW_Pos = BandwidthAnalysis(TN_old_passive_M, 'interp',0,'Pos',.5,'OldM-BW_Pos')
BandwidthAll.Old.M.BW_Neg = BandwidthAnalysis(TN_old_passive_M, 'interp',0,'Neg',.5,'OldM-BW_Neg')


BandwidthAll.Old.F.BW_Pos = BandwidthAnalysis(TN_old_passive_F, 'interp',0,'Pos',.5,'OldF-BW_Pos')
BandwidthAll.Old.F.BW_Neg =BandwidthAnalysis(TN_old_passive_F, 'interp',0,'Neg',.5,'OldF-BW_Neg')



%% ScatterLevels Plots 
% 
% PlotScatterLevels({'Old','Young'},{'Tones','20','10','0'},'OldYoungAll-BW_PosByAnimal' ,BandwidthAll.Old.All.BW_Pos{4},...
%                                  BandwidthAll.Young.All.BW_Pos{4}) 
%  
% PlotScatterLevels({'Old','Young'},{'Tones','20','10','0'},'OldYoungAll-BW_NegByAnimal' ,BandwidthAll.Old.All.BW_Neg{4},...
%                                  BandwidthAll.Young.All.BW_Neg{4}) 
%                                                           
%                              
%  PlotScatterLevels({'Female','Male'},{'Tones','20','10','0'},...
%                                  'YoungMF-BW-ByAnimal_Pos',BandwidthAll.Young.F.BW_Pos{4},...
%                                  BandwidthAll.Young.M.BW_Pos{4});
%                              
%  PlotScatterLevels({'Female','Male'},{'Tones','20','10','0'},...
%                                  'YoungMF-BW-ByAnimal_Neg',BandwidthAll.Young.F.BW_Neg{4},...
%                                  BandwidthAll.Young.M.BW_Neg{4});
%                                          
%                              
%                              
%  PlotScatterLevels( {'Aging Male','Young Male'},{'Tones','20','10','0'},...
%                                  'YoungOldM-BW-ByAnimal_Pos',BandwidthAll.Old.M.BW_Pos{4},...
%                                  BandwidthAll.Young.M.BW_Pos{4});                              
%                              
%  PlotScatterLevels({'Aging Male','Young Male'},{'Tones','20','10','0'},...
%                                  'YoungOldM-BRFS-ByAnimal_Neg',BandwidthAll.Old.M.Neg{4},...
%                                  BandwidthAll.Young.M.Neg{4});     
% 


%% Stats 
  % BRFS

BandwidthAll.Stats.OldYoung.Neg = CompareOldYoung(BandwidthAll.Young.All.Neg,...
                                 BandwidthAll.Old.All.Neg,{'Young','Aging'},...
                                  {'Tones','20','10','0'},...
                                 'OldYoungAll-Neg');
                             
BandwidthAll.Stats.OldYoung.Pos = CompareOldYoung(BandwidthAll.Young.All.Pos,...
                                 BandwidthAll.Old.All.Pos,{'Young','Aging'},......
                                  {'Tones','20','10','0'},...
                                  'OldYoungAll-Pos')

BandwidthAll.Stats.OldYoung.Neg = CompareOldYoung(BandwidthAll.Young.All.Neg,...
                                 BandwidthAll.Old.All.Neg,{'Young','Aging'},...
                                  {'Tones','20','10','0'},...
                                 'OldYoungAll-Neg');
                             
BandwidthAll.Stats.OldYoung.Pos = CompareOldYoung(BandwidthAll.Young.All.Pos,...
                                 BandwidthAll.Old.All.Pos,{'Young','Aging'},......
                                  {'Tones','20','10','0'},...
                                 'OldYoungAll-Pos')






BandwidthAll.Stats.OldYoung.Neg = CompareOldYoung(BandwidthAll.Young.All.Neg,...
                                 BandwidthAll.Old.All.Neg,{'Young','Aging'},...
                                  {'Tones','20','10','0'},...
                                 'OldYoungAll-Neg');
                             
BandwidthAll.Stats.OldYoung.Pos = CompareOldYoung(BandwidthAll.Young.All.Pos,...
                                 BandwidthAll.Old.All.Pos,{'Young','Aging'},......
                                  {'Tones','20','10','0'},...
                                 'OldYoungAll-Pos');   
                         
BandwidthAll.Stats.OldYoung.Neg = CompareOldYoung(BandwidthAll.Young.All.Neg,...
                                 BandwidthAll.Old.All.Neg,{'Young','Aging'},...
                                  {'Tones','20','10','0'},...
                                 'OldYoungAll-Neg');
                             
BandwidthAll.Stats.OldYoung.Pos = CompareOldYoung(BandwidthAll.Young.All.Pos,...
                                 BandwidthAll.Old.All.Pos,{'Young','Aging'},......
                                  {'Tones','20','10','0'},...
                                 'OldYoungAll-Pos'); 
 

                                                      
BandwidthAll.Stats.YoungMF.Pos = CompareOldYoung(BandwidthAll.Young.F.Pos,...
                                 BandwidthAll.Young.M.Pos,{'Female','Male'},...
                                  {'Tones','20','10','0'},...
                                 'YoungMF-Pos');
BandwidthAll.Stats.YoungMF.Neg = CompareOldYoung(BandwidthAll.Young.F.Neg,...
                                 BandwidthAll.Young.M.Neg,{'Female','Male'},...
                                  {'Tones','20','10','0'},...
                                 'YoungMF-Neg');                                

               
BandwidthAll.Stats.YoungOldM.Pos = CompareOldYoung(BandwidthAll.Young.M.Pos,...
                                 BandwidthAll.Old.M.Pos,{'Young-M','Old-M'},...
                                  {'Tones','20','10','0'},...
                                 'YoungOldM-Pos');                               
                             
BandwidthAll.Stats.YoungOldM.Neg = CompareOldYoung(BandwidthAll.Young.M.Neg,...
                                 BandwidthAll.Old.M.Neg,{'Young-M','old-M'},...
                                  {'Tones','20','10','0'},...
                                 'YoungOldM-Neg');     

BandwidthAll.Stats.YoungOldF.Pos = CompareOldYoung(BandwidthAll.Young.F.Pos,...
                                 BandwidthAll.Old.F.Pos,{'Young-F','Old-F'},...
                                  {'Tones','20','10','0'},...
                                 'YoungOldF-Pos');                               
                             
BandwidthAll.Stats.YoungOldF.Neg = CompareOldYoung(BandwidthAll.Young.F.Neg,...
                                 BandwidthAll.Old.F.Neg,{'Young-F','Old-F'},...
                                  {'Tones','20','10','0'},...
                                 'YoungOldF-Neg');     



                                                          
     % Traditional BW

BandwidthAll.Stats.OldYoung.BW_Pos = CompareOldYoung(BandwidthAll.Young.All.BW_Pos,...
                                 BandwidthAll.Old.All.BW_Pos,{'Young','Aging'},......
                                  {'Tones','20','10','0'},...
                                 'OldYoungAll-BW_Pos');
            

BandwidthAll.Stats.OldYoung.BW_Neg = CompareOldYoung(BandwidthAll.Young.All.BW_Neg,...
                                 BandwidthAll.Old.All.BW_Neg,{'Young','Aging'},......
                                  {'Tones','20','10','0'},...
                                 'OldYoungAll-BW_Neg');




BandwidthAll.Stats.OldYoung.BW_Pos = CompareOldYoung(BandwidthAll.Young.All.BW_Pos,...
                                 BandwidthAll.Old.All.BW_Pos,{'Young','Aging'},......
                                  {'Tones','20','10','0'},...
                                 'OldYoungAll-BW_Pos');
            

BandwidthAll.Stats.OldYoung.BW_Neg = CompareOldYoung(BandwidthAll.Young.All.BW_Neg,...
                                 BandwidthAll.Old.All.BW_Neg,{'Young','Aging'},......
                                  {'Tones','20','10','0'},...
                                 'OldYoungAll-BW_Neg');




                             
BandwidthAll.Stats.YoungMF.BW_Pos = CompareOldYoung(BandwidthAll.Young.M.BW_Pos,...
                                 BandwidthAll.Young.F.BW_Pos,{'Male','Female'}, ...
                                  {'Tones','20','10','0'},...
                                 'YoungMF-BW_Pos');
                             
BandwidthAll.Stats.YoungMF.BW_Neg = CompareOldYoung(BandwidthAll.Young.M.BW_Neg,...
                                 BandwidthAll.Young.F.BW_Neg,{'Male','Female'},...
                                  {'Tones','20','10','0'},...
                                 'YoungMF-BW_Neg');            
                             
                             
BandwidthAll.Stats.YoungOldM.BW_Pos = CompareOldYoung(BandwidthAll.Old.M.BW_Pos,...
                                 BandwidthAll.Young.M.BW_Pos,{'Old-M','Young-M'}, ...
                                  {'Tones','20','10','0'},...
                                 'YoungOldM-BW_Pos');                               
                             
BandwidthAll.Stats.YoungOldM.BW_Neg = CompareOldYoung(BandwidthAll.Old.M.BW_Neg,...
                                 BandwidthAll.Young.M.BW_Neg,{'Old-M','Young-M'},...
                                  {'Tones','20','10','0'},...
                                 'YoungOldM-BW_Neg');  
                              
                     
                                                               
  BandwidthAll.Stats.YoungOldF.BW_Pos = CompareOldYoung(BandwidthAll.Young.F.BW_Pos,...
                                 BandwidthAll.Old.F.BW_Pos,{'Young-F','Old-F'},...
                                  {'Tones','20','10','0'},...
                                 'YoungOldF-BW_Pos');  
                             
                              
 BandwidthAll.Stats.YoungOldF.BW_Neg = CompareOldYoung(BandwidthAll.Young.F.BW_Neg,...
                                 BandwidthAll.Old.F.BW_Neg,{'Young-F','Old-F'},...
                                  {'Tones','20','10','0'},...
                                 'YoungOldF-BW_Neg');  
       
                                     
                            
                             
                             
save('BandwidthStats','BandwidthAll')



%% Correlations 
Corr_results  = []

Corr_results.YoungGood = CorrelationsActivePassive(TNBehavior( (~old_idx) & (dPrime_gt1_idx) ,:),'YoungGood_m') ;
Corr_results.YoungBad = CorrelationsActivePassive(TNBehavior( (~old_idx) & (~dPrime_gt1_idx) ,:),'YoungBad_m' );
Corr_results.OldGood = CorrelationsActivePassive(TNBehavior( (old_idx) & (dPrime_gt1_idx) ,:),'OldGood_m');
Corr_results.OldBad = CorrelationsActivePassive(TNBehavior( (old_idx) & (~dPrime_gt1_idx) ,:),'OldBad_m' );

 CorrelationsActivePassive(TNBehavior( (~old_idx) & (dPrime_gt1_idx) ,:),'YoungGood_m','Miss') ;
 CorrelationsActivePassive(TNBehavior( (~old_idx) & (~dPrime_gt1_idx) ,:),'YoungBad_m','Miss' );
 CorrelationsActivePassive(TNBehavior( (old_idx) & (dPrime_gt1_idx) ,:),'OldGood_m','Miss');
 CorrelationsActivePassive(TNBehavior( (old_idx) & (~dPrime_gt1_idx) ,:),'OldBad_m','Miss' );

 CorrelationsActivePassive(TNBehavior( (~old_idx) & (dPrime_gt1_idx) ,:),'YoungGood_e','Early') ;
 CorrelationsActivePassive(TNBehavior( (~old_idx) & (~dPrime_gt1_idx) ,:),'YoungBad_e','Early' );
 CorrelationsActivePassive(TNBehavior( (old_idx) & (dPrime_gt1_idx) ,:),'OldGood_e','Early');
 CorrelationsActivePassive(TNBehavior( (old_idx) & (~dPrime_gt1_idx) ,:),'OldBad_e','Early' );



TN_young_good_p.CorrByAnimal = CorrelationsByAnimal(TN_young_good_p,1);
TN_young_bad_p.CorrByAnimal = CorrelationsByAnimal(TN_young_bad_p,1);
TN_old_good_p.CorrByAnimal = CorrelationsByAnimal(TN_old_good_p,1);
TN_old_bad_p.CorrByAnimal = CorrelationsByAnimal(TN_old_bad_p,1);

 CorrsByDistance(TN_young_good_p,10,'YoungGood_passive')
 CorrsByDistance(TN_young_bad_p,10,'YoungBad_passive')
 CorrsByDistance(TN_old_good_p,10,'OldGood_passive')
 CorrsByDistance(TN_old_bad_p,10,'OldBad_passive')
 
 CorrsByDistance(TNBehavior(dPrime_gt1_idx &~old_idx,2),10,'YoungGood_active')
 CorrsByDistance(TNBehavior((~dPrime_gt1_idx)&~old_idx,2),10,'YoungBad_active')
 CorrsByDistance(TNBehavior( dPrime_gt1_idx & old_idx,2),10,'OldGood_active')
 CorrsByDistance(TNBehavior((~dPrime_gt1_idx)& old_idx,2),10,'OldBad_active')

 TN_young_passive_F.CorrByAnimal = CorrelationsByAnimal(TN_young_passive_F,1);
 TN_young_passive_M.CorrByAnimal = CorrelationsByAnimal(TN_young_passive_M,1);
 TN_old_passive_M.CorrByAnimal = CorrelationsByAnimal(TN_old_passive_M,1);
 TN_old_passive_F.CorrByAnimal =  CorrelationsByAnimal(TN_old_passive_F,1);
 TN_young_passive.CorrByAnimal = CorrelationsByAnimal(TN_young_passive,1); 
 TN_old_passive.CorrByAnimal = CorrelationsByAnimal(TN_old_passive,1);
 %Compare2Anova(g1,g2,GroupNames,levels,SaveName,varnames,full_levels)

CorrDist.YoungAll = CorrsByDistance(TN_young_passive,10,'YoungAll')
CorrDist.YoungM = CorrsByDistance(TN_young_passive_M,10,'YoungMale')
CorrDist.YoungF = CorrsByDistance(TN_young_passive_F,10,'YoungFemale')
CorrDist.AgingAll = CorrsByDistance(TN_old_passive,10,'AgingAll')
CorrDist.AgingM = CorrsByDistance(TN_old_passive_M,10,'AgingMale')
CorrDist.AgingF = CorrsByDistance(TN_old_passive_F,10,'AgingFemale')
 
% Active Correlations 
ActiveCorrs.Young.All= CorrelationsByAnimal(TNBehaviorYoung(:,2),1); 
ActiveCorrs.Young.M = CorrelationsByAnimal(TNBehaviorYoung_Male(:,2),1);
ActiveCorrs.Young.F = CorrelationsByAnimal(TNBehaviorYoung_Female(:,2),1); 
ActiveCorrs.Old.All = CorrelationsByAnimal(TNBehaviorOld(:,2),1); 
ActiveCorrs.Old.M = CorrelationsByAnimal(TNBehaviorOld_Male(:,2),1); 
ActiveCorrs.Old.F = CorrelationsByAnimal(TNBehaviorOld_Female(:,2),1); 

CorrDist.Active.YoungAll = CorrsByDistance(TNBehaviorYoung(:,2),10,'YoungAll')
CorrDist.Active.YoungM = CorrsByDistance(TNBehaviorYoung_Male(:,2),10,'YoungMale')
CorrDist.Active.YoungF = CorrsByDistance(TNBehaviorYoung_Female(:,2),10,'YoungFemale')
CorrDist.Active.AgingAll = CorrsByDistance(TNBehaviorOld(:,2),10,'AgingAll')
CorrDist.Active.AgingM = CorrsByDistance(TNBehaviorOld_Male(:,2),10,'AgingMale')
CorrDist.Active.AgingF = CorrsByDistance(TNBehaviorOld_Female(:,2),10,'AgingFemale')


CompareCorrsByDistance(CorrDist.YoungAll,CorrDist.AgingAll,'YoungAgingAll')
CompareCorrsByDistance(CorrDist.YoungM,CorrDist.AgingM,'YoungAgingMale')
CompareCorrsByDistance(CorrDist.YoungF,CorrDist.AgingF,'YoungAgingFemale')

 
 
 Corr.YoungMF.SignalCorr = CorrelationsAnalysisByCell(TN_young_passive_M.CorrByAnimal.LCorr',...
TN_young_passive_F.CorrByAnimal.LCorr',{'Male','Female'},{'Tones','20','10','0'},...
'SignalCorr-Cell') ;
                              
 Corr.YoungMF.NoiseCorr =  CorrelationsAnalysisByCell(TN_young_passive_M.CorrByAnimal.NCorr',...
TN_young_passive_F.CorrByAnimal.NCorr',{'Male','Female'},{'Tones','20','10','0'},...
'NoiseCorr-Cell') ;




 Corr.YoungOldM.SignalCorr =  CorrelationsAnalysisByCell(TN_young_passive_M.CorrByAnimal.LCorr',...
TN_old_passive_M.CorrByAnimal.LCorr',{'Males-Young','Old'},{'Tones','20','10','0'},...
'SignalCorr-Cell') ;

Corr.YoungOldM.NoiseCorr =  CorrelationsAnalysisByCell(TN_young_passive_M.CorrByAnimal.NCorr',...
TN_old_passive_M.CorrByAnimal.NCorr',{'Males-Young','Old'},{'Tones','20','10','0'},...
'NoiseCorr-Cell') ;

Corr.YoungOldF.SignalCorr =  CorrelationsAnalysisByCell(TN_young_passive_F.CorrByAnimal.LCorr',...
TN_old_passive_F.CorrByAnimal.LCorr',{'Females-Young','Old'},{'Tones','20','10','0'},...
'SignalCorr-Cell') ;

Corr.YoungOldF.NoiseCorr =  CorrelationsAnalysisByCell(TN_young_passive_F.CorrByAnimal.NCorr',...
TN_old_passive_F.CorrByAnimal.NCorr',{'Females-Young','Old'},{'Tones','20','10','0'},...
'NoiseCorr-Cell') ;


% Corr.OldYoung.SignalCorr = CorrelationsAnalysisByCell(TN_young_passive2.CorrByAnimal.LCorr',...
% TN_old_passive2.CorrByAnimal.LCorr',{'Young','Aging'},{'Tones','20','10','0'},...
% 'SignalCorr-Cell') ;
% 
% Corr.OldYoung.NoiseCorr =  CorrelationsAnalysisByCell(TN_young_passive2.CorrByAnimal.NCorr',...
% TN_old_passive2.CorrByAnimal.NCorr',{'Young','Aging'},{'Tones','20','10','0'},...
% 'NoiseCorr-Cell') ;



Corr.OldYoung.SignalCorr = CorrelationsAnalysisByCell(TN_young_passive2.CorrByAnimal.LCorr',...
TN_old_passive.CorrByAnimal.LCorr',{'Young','Aging'},{'Tones','20','10','0'},...
'SignalCorr-Cell') ;

Corr.OldYoung.NoiseCorr =  CorrelationsAnalysisByCell(TN_young_passive2.CorrByAnimal.NCorr',...
TN_old_passive.CorrByAnimal.NCorr',{'Young','Aging'},{'Tones','20','10','0'},...
'NoiseCorr-Cell') ;


 
  Corr.YoungMF.SignalCorrAnimal = CorrelationsAnalysisByAnimal(TN_young_passive_M.CorrByAnimal.LCorrAnimal',...
                                  TN_young_passive_F.CorrByAnimal.LCorrAnimal',...
                                  'SignalCorr-byAnimal',{'Male','Female'},...
                                  {'Tones','20','10','0'}) ;  
                              
 Corr.YoungMF.NoiseCorrAnimal = CorrelationsAnalysisByAnimal(TN_young_passive_M.CorrByAnimal.NCorrAnimal',...
                                  TN_young_passive_F.CorrByAnimal.NCorrAnimal',...
                                  'NoiseCorr-byAnimal',{'Male','Female'},...
                                  {'Tones','20','10','0'}) ;                           
                                                                                                                                      
 
                           
  Corr.YoungOldM.SignalCorrAnimal = CorrelationsAnalysisByAnimal(TN_young_passive_M.CorrByAnimal.LCorrAnimal',...
                                  TN_old_passive_M.CorrByAnimal.LCorrAnimal',...
                                  'SignalCorr-byAnimal',{'Males-Young','Aging'},...
                                  {'Tones','20','10','0'}) ;                           
                                                     
  
   Corr.YoungOldM.NoiseCorrAnimal = CorrelationsAnalysisByAnimal(TN_young_passive_M.CorrByAnimal.NCorrAnimal',...
                                  TN_old_passive_M.CorrByAnimal.NCorrAnimal',...
                                  'NoiseCorr-byAnimal',{'Males-Young','Aging'},...
                                  {'Tones','20','10','0'}) ;
                              
                              
    Corr.YoungOldF.SignalCorrAnimal = CorrelationsAnalysisByAnimal(TN_young_passive_F.CorrByAnimal.LCorrAnimal',...
                                  TN_old_passive_F.CorrByAnimal.LCorrAnimal',...
                                  'SignalCorr-byAnimal',{'Females-Young','Aging'},...
                                  {'Tones','20','10','0'}) ;                           
                                                     
  
   Corr.YoungOldF.NoiseCorrAnimal = CorrelationsAnalysisByAnimal(TN_young_passive_F.CorrByAnimal.NCorrAnimal',...
                                  TN_old_passive_F.CorrByAnimal.NCorrAnimal',...
                                  'NoiseCorr-byAnimal',{'Females-Young','Aging'},...
                                  {'Tones','20','10','0'}) ;                             
                                                                                                                                           
                              
                              
  Corr.OldYoung.SignalCorrAnimal = CorrelationsAnalysisByAnimal(TN_young_passive2.CorrByAnimal.LCorrAnimal',...
                                  TN_old_passive2.CorrByAnimal.LCorrAnimal',...
                                  'SignalCorr-byAnimal',{'Young','Aging'},...
                                  {'Tones','20','10','0'}) ;                           
  Corr.OldYoung.NoiseCorrAnimal = CorrelationsAnalysisByAnimal(TN_young_passive2.CorrByAnimal.NCorrAnimal',...
                                  TN_old_passive2.CorrByAnimal.NCorrAnimal',...
                                  'NoiseCorr-byAnimal',{'Young','Aging'},...
                                  {'Tones','20','10','0'}) ; 
  
  Corr.OldYoung.SignalCorrAnimal = CorrelationsAnalysisByAnimal(TN_young_passive2.CorrByAnimal.LCorrAnimal',...
                                  TN_old_passive.CorrByAnimal.LCorrAnimal',...
                                  'SignalCorr-byAnimal',{'Young','Aging'},...
                                  {'Tones','20','10','0'}) ;                           
  Corr.OldYoung.NoiseCorrAnimal = CorrelationsAnalysisByAnimal(TN_young_passive2.CorrByAnimal.NCorrAnimal',...
                                  TN_old_passive.CorrByAnimal.NCorrAnimal',...
                                  'NoiseCorr-byAnimal',{'Young','Aging'},...
                                  {'Tones','20','10','0'}) ;   
 
                        



                     

     
                              

  % Plots that remove the effect of Group to focus on the effects of level
 Corr.OldYoung.NoiseDiffCell = CorrelationsAnalysisByAnimal(TN_young_passive2.CorrByAnimal.NCorrDiff',...
                                  TN_old_passive.CorrByAnimal.NCorrDiff',...
                                  'NoiseDiff-byCell',{'Young','Aging-Rev'},...
                                  {'Tones','20','10','0'}) ;                           
                                                    
   Corr.OldYoung.SignalDiffCell = CorrelationsAnalysisByAnimal(TN_young_passive2.CorrByAnimal.LCorrDiff',...
                                  TN_old_passive.CorrByAnimal.LCorrDiff',...
                                  'SignalDiff-byCell',{'Young','Aging-Rev'},...
                                  {'Tones','20','10','0'}) ;                           
                                                       

 Corr.OldYoungM.NoiseDiffCell = CorrelationsAnalysisByAnimal(TN_young_passive_M.CorrByAnimal.NCorrDiff',...
                                  TN_old_passive_M.CorrByAnimal.NCorrDiff',...
                                  'NoiseDiff-byCell',{'Young','Aging-Rev'},...
                                  {'Tones','20','10','0'}) ;                           
                                                    
   Corr.OldYoungM.SignalDiffCell = CorrelationsAnalysisByAnimal(TN_young_passive2.CorrByAnimal.LCorrDiff',...
                                  TN_old_passive_M.CorrByAnimal.LCorrDiff',...
                                  'SignalDiff-byCell',{'YoungM','AgingM'},...
                                  {'Tones','20','10','0'}) ;                           
                                                       





                              
  Corr.OldYoung.NoiseDiffAnimal = CorrelationsAnalysisByAnimal(TN_young_passive2.CorrByAnimal.NCorrDiffAnimal',...
                                  TN_old_passive2.CorrByAnimal.NCorrDiffAnimal',...
                                  'NoiseDiff-byAnimal',{'Young','Aging-Rev'},...
                                  {'Tones','20','10','0'}) ;                           
                                                    
   Corr.OldYoung.SignalDiffAnimal = CorrelationsAnalysisByAnimal(TN_young_passive2.CorrByAnimal.LCorrDiffAnimal',...
                                  TN_old_passive2.CorrByAnimal.LCorrDiffAnimal',...
                                  'SignalDiff-byAnimal',{'Young','Aging-Rev'},...
                                  {'Tones','20','10','0'}) ;                           
                                                    
 
             
  Corr.OldYoungM.NoiseDiffAnimal = CorrelationsAnalysisByAnimal(TN_young_passive_M.CorrByAnimal.NCorrDiffAnimal',...
                                  TN_old_passive_M.CorrByAnimal.NCorrDiffAnimal',...
                                  'NoiseDiff-byAnimal',{'YoungM','AgingM'},...
                                  {'Tones','20','10','0'}) ;                           
                                                    
   Corr.OldYoungM.SignalDiffAnimal = CorrelationsAnalysisByAnimal(TN_young_passive_M.CorrByAnimal.LCorrDiffAnimal',...
                                  TN_old_passive_M.CorrByAnimal.LCorrDiffAnimal',...
                                  'SignalDiff-byAnimal',{'YoungM','AgingM'},...
                                  {'Tones','20','10','0'}) ;                           
                              
 
    % Correlation by distance                    
                              
                    


%% Cluster comparision
All_Passive = Fluoro_to_Table(TNBehavior(:,3) ,49000);

All_Passive.Clusters =Cluster_DF(All_Passive,'K-means',20);

TN_old_good_p.Combined_Classes =   getClusters(All_Passive,old_idx & dPrime_gt1_idx);
TN_old_bad_p.Combined_Classes =   getClusters(All_Passive,old_idx & (~dPrime_gt1_idx) );
TN_young_good_p.Combined_Classes = getClusters(All_Passive, (~old_idx) & (dPrime_gt1_idx));
TN_young_bad_p.Combined_Classes =  getClusters(All_Passive, (~old_idx) & (~dPrime_gt1_idx ));


Plot_Clusters(All_Passive,All_Passive.Clusters.Clusters);

Plot_Clusters(TN_young_good_p,TN_young_good_p.Combined_Classes,'Young_Good')
Plot_Clusters(TN_young_bad_p,TN_young_bad_p.Combined_Classes,'Young_Bad')
Plot_Clusters(TN_old_good_p,TN_old_good_p.Combined_Classes,'Aging_Good')
Plot_Clusters(TN_old_bad_p,TN_old_bad_p.Combined_Classes,'Aging_Bad')


Diversity = struct()


Diversity.Stats.Young_GB = Diversityanalysis(TN_young_good_p,TN_young_bad_p,{'YoungGood','YoungBad'});

Diversity.Stats.Good_YO =Diversityanalysis(TN_young_good_p,TN_old_good_p,{'YoungGood','OldGood'});

Diversity.Stats.Bad_YO = Diversityanalysis(TN_young_bad_p,TN_old_bad_p,{'YoungBad','OldBad'});

Diversity.Stats.Old_GB = Diversityanalysis(TN_old_good_p,TN_old_bad_p,{'OldGood','OldBad'});



%%  Bayes 

BayesModels = struct()

 BayesModels.YoungGood = BayesClassiferActive(TNBehavior( (~old_idx) & (dPrime_gt1_idx) ,2));
 BayesModels.YoungBad = BayesClassiferActive(TNBehavior( (~old_idx) & (~dPrime_gt1_idx) ,2));
 BayesModels.OldBad = BayesClassiferActive(TNBehavior( (old_idx) & (~dPrime_gt1_idx) ,2));
 BayesModels.OldGood = BayesClassiferActive(TNBehavior( (old_idx) & (dPrime_gt1_idx) ,2));

BayesModels.YoungGood.Stats = PlotBayesActive(BayesModels.YoungGood,'YoungGood');
BayesModels.YoungBad.Stats = PlotBayesActive(BayesModels.YoungBad,'YoungBad');
BayesModels.OldGood.Stats = PlotBayesActive(BayesModels.OldGood,'OldGood');
BayesModels.OldBad.Stats = PlotBayesActive(BayesModels.OldBad,'OldBad');


PlotBayesFig(BayesModels,'TimeLossTotal_Hits','BayesTiming-Hits')
PlotBayesFig(BayesModels,'NumbersLossTotal_Hits','BayesNumbers-Hits')

PlotBayesFig(BayesModels,'NumbersLossTotal_Hits','BayesNumbers-Hits')
PlotBayesFig(BayesModels,'NumbersLossTotal_HitsEarly','BayesNumbers-HitsEarly')


Compare2AnovaBehavior({BayesModels.YoungGood.Stats.TimeLossTotal_HitsEarly.MaxDeltaPrediction,...
                       BayesModels.YoungBad.Stats.TimeLossTotal_HitsEarly.MaxDeltaPrediction},...
                      {BayesModels.OldGood.Stats.TimeLossTotal_HitsEarly.MaxDeltaPrediction,...
                       BayesModels.OldBad.Stats.TimeLossTotal_HitsEarly.MaxDeltaPrediction},...
                                  {'Young','Aging'},...
                                  {'GoodPerformance','PoorPerformance'},...
                                  'BayesModels.Timing_HitsEarly') ;
                              
                              
Compare2AnovaBehavior({BayesModels.YoungGood.Stats.TimeLossTotal_HitsEarly.MaxPrediction,...
                       BayesModels.YoungBad.Stats.TimeLossTotal_HitsEarly.MaxPrediction},...
                      {BayesModels.OldGood.Stats.TimeLossTotal_HitsEarly.MaxPrediction,...
                       BayesModels.OldBad.Stats.TimeLossTotal_HitsEarly.MaxPrediction},...
                                  {'Young','Aging'},...
                                  {'GoodPerformance','PoorPerformance'},...
                                  'BayesTiming_HitsEarly_Max') ;                               

Compare2AnovaBehavior({BayesModels.YoungGood.Stats.TimeLossTotal_Hits.MaxPrediction,...
                       BayesModels.YoungBad.Stats.TimeLossTotal_Hits.MaxPrediction},...
                      {BayesModels.OldGood.Stats.TimeLossTotal_Hits.MaxPrediction,...
                       BayesModels.OldBad.Stats.TimeLossTotal_Hits.MaxPrediction},...
                                  {'Young','Aging'},...
                                  {'GoodPerformance','PoorPerformance'},...
                                  'BayesTiming_Max') ;                               
 
Compare2AnovaBehavior({BayesModels.YoungGood.Stats.TimeLossTotal_Hits.MaxDeltaPrediction,...
                       BayesModels.YoungBad.Stats.TimeLossTotal_Hits.MaxDeltaPrediction},...
                      {BayesModels.OldGood.Stats.TimeLossTotal_Hits.MaxDeltaPrediction,...
                       BayesModels.OldBad.Stats.TimeLossTotal_Hits.MaxDeltaPrediction},...
                                  {'Young','Aging'},...
                                  {'GoodPerformance','PoorPerformance'},...
                                  'BayesTiming') ;                              
                              
