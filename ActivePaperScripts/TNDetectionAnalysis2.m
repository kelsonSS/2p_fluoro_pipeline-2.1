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
TNBehaviorInfo = getAnimalInfo(TNBehavior)

old_idx = contains(lower(TNBehavior(:,3)),'ia');

TNBehaviorYoung = TNBehavior(~old_idx,:);
TNBehaviorAging = TNBehavior(old_idx,:);
TNAgingAnimalInfo = getAnimalInfo(TNBehaviorAging(:,4));
TNYoungAnimalInfo = getAnimalInfo(TNBehaviorYoung(:,4));


%% create output struct
BehaviorResults = struct();

BehaviorResults.Young.CellInfo = CreateBehaviorCellList(TNBehaviorYoung);
BehaviorResults.Aging.CellInfo = CreateBehaviorCellList(TNBehaviorAging);

AnimalIDS = cellfun(@(x) strsplit(x,filesep),TNBehavior(:,3),'UniformOutput',0);
 AnimalIDS = unique( cellfun(@(x) x{4},AnimalIDS,'UniformOutput',0));
 BehaviorResults.AnimalID = AnimalIDS;
 clear AnimalIDS
 
%% clustering 
% BehaviorResults.Young.Clust.Passive = Cluster_DF(TNBehaviorYoung(:,1),'K-means');
% BehaviorResults.Young.Clust.Active = Cluster_DF(TNBehaviorYoung(:,2),'K-means',9,...
% 'normalized',BehaviorResults.Young.Clust.Passive.Centroids);
% 
% BehaviorResults.Aging.Clust.Passive = Cluster_DF(TNBehaviorAging(:,1),'K-means');
% BehaviorResults.Aging.Clust.Active = Cluster_DF(TNBehaviorAging(:,2),'K-means',9,...
% 'normalized',BehaviorResults.Aging.Clust.Passive.Centroids);
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
%  Aging_passive_clusters = BehaviorResults.Aging.Clust.Passive.Clusters(...
%      BehaviorResults.Aging.Clust.Passive.Neuron);
%  Aging_active_clusters = BehaviorResults.Aging.Clust.Active.Clusters(...
%      BehaviorResults.Aging.Clust.Active.Neuron);
% 
%  % plot Young
%  trans_mat_young = histcounts2(Young_passive_clusters,Young_active_clusters(1:2892));
% trans_mat_young = trans_mat_young ./ sum(trans_mat_young(:)) ;
%  figure; imagesc(trans_mat_young)
 % Plot Aging 
 
 %%
 
  behavior_path ='Z:\Kelson\TNDetectionAnalysis';
% load in behavior 
for ii = 1:length(BehaviorResults.AnimalID)
    curr_id = BehaviorResults.AnimalID{ii};
    
    % if animal is not in path, add it to path
 %   if  exist(fullfile(behavior_path,[curr_id '.mat']),'file')
           BehavioralAnalysisByAnimal(curr_id,'Totals');
           pause
 %   end 
        
      
    if contains(curr_id,'IA','IgnoreCase',true)
        
        BehaviorResults.Aging.Behavior.( sprintf('animal_%s',curr_id) ) = ...
        load(fullfile(behavior_path,curr_id));
    
              
    else 
         BehaviorResults.Young.Behavior.( sprintf('animal_%s',curr_id) ) = ...
        load(fullfile(behavior_path,curr_id));
    end 
    
        
end 

TN_passive =Munge_DF(TNBehasvior(:,1));
TN_active = Munge_DF(TNBehavior(:,2));

% analyze behavior
TNBehaviorResultsYoung
TNBehaviorResultsAging

BehaviorResults.Aging.Behavior.group = MungeBehaviorGroupData(BehaviorResults.Aging.Behavior);
BehaviorResults.Young.Behavior.group = MungeBehaviorGroupData(BehaviorResults.Young.Behavior);

%Plot group Behavior 
PlotGroupedData(BehaviorResults.Aging.Behavior.group)
PlotGroupedData(BehaviorResults.Young.Behavior.group)


% analyse Passive Data
TN_young_passive = Fluoro_to_Table(TNBehaviorYoung(:,3),50000);
TN_young_passive.AnimalInfo = getAnimalInfo(TN_young_passive.DataDirs);
TN_old_passive =  Fluoro_to_Table(TNBehaviorAging(:,3),50000);
TN_old_passive.AnimalInfo = getAnimalInfo(TN_old_passive.DataDirs);


 % generate datasets seperated by sex
  old_F_idx = contains(TN_old_passive.AnimalInfo{:,2},'F'); 
  young_F_idx= contains(TN_young_passive.AnimalInfo{:,2},'F');
 
  
  TN_old_passive_F = Fluoro_to_Table(TN_old_passive.DataDirs(old_F_idx),50000);
   TN_old_passive_M = Fluoro_to_Table(TN_old_passive.DataDirs(~old_F_idx),50000);
% 
%   
    TN_young_passive_F = Fluoro_to_Table(TN_young_passive.DataDirs(young_F_idx),50000);
   TN_young_passive_M = Fluoro_to_Table(TN_young_passive.DataDirs(~young_F_idx),50000);
%   


%Intensity
 %% Intensity
   Intensity = struct()
  % Age Intensity

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
  
  
   Intensity.Stats.AgingYoungM.Mean = FluoroIntensityAnalysis(TN_young_passive_M,TN_old_passive_M,...
      'mean',{'Young-M','Aging-M'})
  
   Intensity.Stats.AgingYoungM.Max =  FluoroIntensityAnalysis(TN_young_passive_M,TN_old_passive_M,...
      'max',{'Young-M','Aging-M'})
  

    % old F vs young F
  
  
   Intensity.Stats.AgingYoungF.Mean = FluoroIntensityAnalysis(TN_young_passive_F,TN_old_passive_F,...
      'mean',{'Young-F','Aging-F'})
  
   Intensity.Stats.AgingYoungF.Max =  FluoroIntensityAnalysis(TN_young_passive_F,TN_old_passive_F,...
      'max',{'Young-F','Aging-F'})
  
  
  
save 'IntensityStats.mat' 'Intensity' 

 %% Timing
    %% Timing Noise 
   
Timings = struct();
Timings.Passive.YoungF = TemporalAnalysisByAnimal(TN_young_passive_F,'all','Timing_young_F_noise');
    
Timings.Passive.YoungM = TemporalAnalysisByAnimal(TN_young_passive_M,'all','Timing_young_M_noise');

Timings.Passive.AgingF = TemporalAnalysisByAnimal(TN_old_passive_F,'all','Timing_Aging_F_noise');

Timings.Passive.AgingM =TemporalAnalysisByAnimal(TN_old_passive_M,'all','Timing_Aging_M_noise');


Timings.Passive.Aging = TemporalAnalysisByAnimal(TN_old_passive,'all','Timing_Aging_All_noise');

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
%       % Aging_M vs Young_M
%       Compareare2Anova(young_timing_noise_M',old_timing_noise_M',{'Males-Young','Aging'},[1:size(old_timing_noise_M,2])
%       
% % Aging_F vs Young_F
%       CompareAgingYoung(young_timing_noise_F',old_timing_noise_F',[1:size(old_timing_noise_M,2)])
%       
      


save('TimingStats','Timings')


%% Bandwidth Comparision
BandwithAll = struct()



% Using BRFS 
BandwidthAll.Young.All.Pos = BandwidthAnalysis(TN_young_passive, 'BRFS',0,'Pos',.5,'YoungAll-Pos')
BandwidthAll.Young.All.Neg = BandwidthAnalysis(TN_young_passive, 'BRFS',0,'Neg',.5,'YoungAll-Neg')

BandwidthAll.Aging.All.Pos = BandwidthAnalysis(TN_old_passive, 'BRFS',0,'Pos',.5,'AgingAll-Pos')
BandwidthAll.Aging.All.Neg =BandwidthAnalysis(TN_old_passive, 'BRFS',0,'Neg',.5,'AgingAll-Neg')

BandwidthAll.Aging.All.Pos = BandwidthAnalysis(TN_old_passive, 'BRFS',0,'Pos',.5,'AgingAll-Pos')
BandwidthAll.Aging.All.Neg =BandwidthAnalysis(TN_old_passive, 'BRFS',0,'Neg',.5,'AgingAll-Neg')



BandwidthAll.Young.F.Pos = BandwidthAnalysis(TN_young_passive_F, 'BRFS',0,'Pos',.5,'YoungF-Pos')
BandwidthAll.Young.F.Neg =BandwidthAnalysis(TN_young_passive_F, 'BRFS',0,'Neg',.5,'YoungF-Neg')

BandwidthAll.Young.M.Pos = BandwidthAnalysis(TN_young_passive_M, 'BRFS',0,'Pos',.5,'YoungM-Pos')
BandwidthAll.Young.M.Neg =BandwidthAnalysis(TN_young_passive_M, 'BRFS',0,'Neg',.5,'YoungM-Neg')


BandwidthAll.Aging.F.Pos = BandwidthAnalysis(TN_old_passive_F, 'BRFS',0,'Pos',.5,'AgingF-Pos')
BandwidthAll.Aging.F.Neg =BandwidthAnalysis(TN_old_passive_F, 'BRFS',0,'Neg',.5,'AgingF-Neg')

BandwidthAll.Aging.M.Pos = BandwidthAnalysis(TN_old_passive_M, 'BRFS',0,'Pos',.5,'AgingM-Pos')
BandwidthAll.Aging.M.Neg =BandwidthAnalysis(TN_old_passive_M, 'BRFS',0,'Neg',.5,'AgingM-Neg')


% using traditional_BW
BandwidthAll.Young.All.BW_Pos = BandwidthAnalysis(TN_young_passive, 'interp',0,'Pos',.5,'YoungAll-BW-Pos')
BandwidthAll.Young.All.BW_Neg =BandwidthAnalysis(TN_young_passive, 'interp',0,'Neg',.5,'YoungAll-BW-Neg')

BandwidthAll.Aging.All.BW_Pos = BandwidthAnalysis(TN_old_passive, 'interp',0,'Pos',.5,'AgingAll-BW-Pos')
BandwidthAll.Aging.All.BW_Neg =BandwidthAnalysis(TN_old_passive, 'interp',0,'Neg',.5,'AgingAll-BW-Neg')

BandwidthAll.Aging.All.BW_Pos = BandwidthAnalysis(TN_old_passive, 'interp',0,'Pos',.5,'AgingAll-BW-Pos')
BandwidthAll.Aging.All.BW_Neg =BandwidthAnalysis(TN_old_passive, 'interp',0,'Neg',.5,'AgingAll-BW-Neg')



BandwidthAll.Young.F.BW_Pos = BandwidthAnalysis(TN_young_passive_F, 'interp',0,'Pos',.5,'YoungF-BW_Pos')
BandwidthAll.Young.F.BW_Neg =BandwidthAnalysis(TN_young_passive_F, 'interp',0,'Neg',.5,'YoungF-BW_Neg')

BandwidthAll.Young.M.BW_Pos = BandwidthAnalysis(TN_young_passive_M, 'interp',0,'Pos',.5,'YoungM-BW_Pos')
BandwidthAll.Young.M.BW_Neg =BandwidthAnalysis(TN_young_passive_M, 'interp',0,'Neg',.5,'YoungM-BW_Neg')

BandwidthAll.Aging.M.BW_Pos = BandwidthAnalysis(TN_old_passive_M, 'interp',0,'Pos',.5,'AgingM-BW_Pos')
BandwidthAll.Aging.M.BW_Neg = BandwidthAnalysis(TN_old_passive_M, 'interp',0,'Neg',.5,'AgingM-BW_Neg')


BandwidthAll.Aging.F.BW_Pos = BandwidthAnalysis(TN_old_passive_F, 'interp',0,'Pos',.5,'AgingF-BW_Pos')
BandwidthAll.Aging.F.BW_Neg =BandwidthAnalysis(TN_old_passive_F, 'interp',0,'Neg',.5,'AgingF-BW_Neg')



%% ScatterLevels Plots 
% 
% PlotScatterLevels({'Aging','Young'},{'Tones','20','10','0'},'AgingYoungAll-BW_PosByAnimal' ,BandwidthAll.Aging.All.BW_Pos{4},...
%                                  BandwidthAll.Young.All.BW_Pos{4}) 
%  
% PlotScatterLevels({'Aging','Young'},{'Tones','20','10','0'},'AgingYoungAll-BW_NegByAnimal' ,BandwidthAll.Aging.All.BW_Neg{4},...
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
%                                  'YoungAgingM-BW-ByAnimal_Pos',BandwidthAll.Aging.M.BW_Pos{4},...
%                                  BandwidthAll.Young.M.BW_Pos{4});                              
%                              
%  PlotScatterLevels({'Aging Male','Young Male'},{'Tones','20','10','0'},...
%                                  'YoungAgingM-BRFS-ByAnimal_Neg',BandwidthAll.Aging.M.Neg{4},...
%                                  BandwidthAll.Young.M.Neg{4});     
% 


%% Stats 
  % BRFS

BandwidthAll.Stats.AgingYoung.Neg = CompareAgingYoung(BandwidthAll.Young.All.Neg,...
                                 BandwidthAll.Aging.All.Neg,{'Young','Aging'},...
                                  {'Tones','20','10','0'},...
                                 'AgingYoungAll-Neg');
                             
BandwidthAll.Stats.AgingYoung.Pos = CompareAgingYoung(BandwidthAll.Young.All.Pos,...
                                 BandwidthAll.Aging.All.Pos,{'Young','Aging'},......
                                  {'Tones','20','10','0'},...
                                  'AgingYoungAll-Pos')

BandwidthAll.Stats.AgingYoung.Neg = CompareAgingYoung(BandwidthAll.Young.All.Neg,...
                                 BandwidthAll.Aging.All.Neg,{'Young','Aging'},...
                                  {'Tones','20','10','0'},...
                                 'AgingYoungAll-Neg');
                             
BandwidthAll.Stats.AgingYoung.Pos = CompareAgingYoung(BandwidthAll.Young.All.Pos,...
                                 BandwidthAll.Aging.All.Pos,{'Young','Aging'},......
                                  {'Tones','20','10','0'},...
                                 'AgingYoungAll-Pos')






BandwidthAll.Stats.AgingYoung.Neg = CompareAgingYoung(BandwidthAll.Young.All.Neg,...
                                 BandwidthAll.Aging.All.Neg,{'Young','Aging'},...
                                  {'Tones','20','10','0'},...
                                 'AgingYoungAll-Neg');
                             
BandwidthAll.Stats.AgingYoung.Pos = CompareAgingYoung(BandwidthAll.Young.All.Pos,...
                                 BandwidthAll.Aging.All.Pos,{'Young','Aging'},......
                                  {'Tones','20','10','0'},...
                                 'AgingYoungAll-Pos');   
                         
BandwidthAll.Stats.AgingYoung.Neg = CompareAgingYoung(BandwidthAll.Young.All.Neg,...
                                 BandwidthAll.Aging.All.Neg,{'Young','Aging'},...
                                  {'Tones','20','10','0'},...
                                 'AgingYoungAll-Neg');
                             
BandwidthAll.Stats.AgingYoung.Pos = CompareAgingYoung(BandwidthAll.Young.All.Pos,...
                                 BandwidthAll.Aging.All.Pos,{'Young','Aging'},......
                                  {'Tones','20','10','0'},...
                                 'AgingYoungAll-Pos'); 
 

                                                      
BandwidthAll.Stats.YoungMF.Pos = CompareAgingYoung(BandwidthAll.Young.F.Pos,...
                                 BandwidthAll.Young.M.Pos,{'Female','Male'},...
                                  {'Tones','20','10','0'},...
                                 'YoungMF-Pos');
BandwidthAll.Stats.YoungMF.Neg = CompareAgingYoung(BandwidthAll.Young.F.Neg,...
                                 BandwidthAll.Young.M.Neg,{'Female','Male'},...
                                  {'Tones','20','10','0'},...
                                 'YoungMF-Neg');                                

               
BandwidthAll.Stats.YoungAgingM.Pos = CompareAgingYoung(BandwidthAll.Young.M.Pos,...
                                 BandwidthAll.Aging.M.Pos,{'Young-M','Aging-M'},...
                                  {'Tones','20','10','0'},...
                                 'YoungAgingM-Pos');                               
                             
BandwidthAll.Stats.YoungAgingM.Neg = CompareAgingYoung(BandwidthAll.Young.M.Neg,...
                                 BandwidthAll.Aging.M.Neg,{'Young-M','old-M'},...
                                  {'Tones','20','10','0'},...
                                 'YoungAgingM-Neg');     

BandwidthAll.Stats.YoungAgingF.Pos = CompareAgingYoung(BandwidthAll.Young.F.Pos,...
                                 BandwidthAll.Aging.F.Pos,{'Young-F','Aging-F'},...
                                  {'Tones','20','10','0'},...
                                 'YoungAgingF-Pos');                               
                             
BandwidthAll.Stats.YoungAgingF.Neg = CompareAgingYoung(BandwidthAll.Young.F.Neg,...
                                 BandwidthAll.Aging.F.Neg,{'Young-F','Aging-F'},...
                                  {'Tones','20','10','0'},...
                                 'YoungAgingF-Neg');     



                                                          
     % Traditional BW

BandwidthAll.Stats.AgingYoung.BW_Pos = CompareAgingYoung(BandwidthAll.Young.All.BW_Pos,...
                                 BandwidthAll.Aging.All.BW_Pos,{'Young','Aging'},......
                                  {'Tones','20','10','0'},...
                                 'AgingYoungAll-BW_Pos');
            

BandwidthAll.Stats.AgingYoung.BW_Neg = CompareAgingYoung(BandwidthAll.Young.All.BW_Neg,...
                                 BandwidthAll.Aging.All.BW_Neg,{'Young','Aging'},......
                                  {'Tones','20','10','0'},...
                                 'AgingYoungAll-BW_Neg');




BandwidthAll.Stats.AgingYoung.BW_Pos = CompareAgingYoung(BandwidthAll.Young.All.BW_Pos,...
                                 BandwidthAll.Aging.All.BW_Pos,{'Young','Aging'},......
                                  {'Tones','20','10','0'},...
                                 'AgingYoungAll-BW_Pos');
            

BandwidthAll.Stats.AgingYoung.BW_Neg = CompareAgingYoung(BandwidthAll.Young.All.BW_Neg,...
                                 BandwidthAll.Aging.All.BW_Neg,{'Young','Aging'},......
                                  {'Tones','20','10','0'},...
                                 'AgingYoungAll-BW_Neg');




                             
BandwidthAll.Stats.YoungMF.BW_Pos = CompareAgingYoung(BandwidthAll.Young.M.BW_Pos,...
                                 BandwidthAll.Young.F.BW_Pos,{'Male','Female'}, ...
                                  {'Tones','20','10','0'},...
                                 'YoungMF-BW_Pos');
                             
BandwidthAll.Stats.YoungMF.BW_Neg = CompareAgingYoung(BandwidthAll.Young.M.BW_Neg,...
                                 BandwidthAll.Young.F.BW_Neg,{'Male','Female'},...
                                  {'Tones','20','10','0'},...
                                 'YoungMF-BW_Neg');            
                             
                             
BandwidthAll.Stats.YoungAgingM.BW_Pos = CompareAgingYoung(BandwidthAll.Aging.M.BW_Pos,...
                                 BandwidthAll.Young.M.BW_Pos,{'Aging-M','Young-M'}, ...
                                  {'Tones','20','10','0'},...
                                 'YoungAgingM-BW_Pos');                               
                             
BandwidthAll.Stats.YoungAgingM.BW_Neg = CompareAgingYoung(BandwidthAll.Aging.M.BW_Neg,...
                                 BandwidthAll.Young.M.BW_Neg,{'Aging-M','Young-M'},...
                                  {'Tones','20','10','0'},...
                                 'YoungAgingM-BW_Neg');  
                              
                     
                                                               
  BandwidthAll.Stats.YoungAgingF.BW_Pos = CompareAgingYoung(BandwidthAll.Young.F.BW_Pos,...
                                 BandwidthAll.Aging.F.BW_Pos,{'Young-F','Aging-F'},...
                                  {'Tones','20','10','0'},...
                                 'YoungAgingF-BW_Pos');  
                             
                              
 BandwidthAll.Stats.YoungAgingF.BW_Neg = CompareAgingYoung(BandwidthAll.Young.F.BW_Neg,...
                                 BandwidthAll.Aging.F.BW_Neg,{'Young-F','Aging-F'},...
                                  {'Tones','20','10','0'},...
                                 'YoungAgingF-BW_Neg');  
       
                                     
                            
                             
                             
save('BandwidthStats','BandwidthAll')



%% Correlations 

 TN_young_passive_F.CorrByAnimal = CorrelationsByAnimal(TN_young_passive_F,1);
 TN_young_passive_M.CorrByAnimal = CorrelationsByAnimal(TN_young_passive_M,1);
 TN_old_passive_M.CorrByAnimal = CorrelationsByAnimal(TN_old_passive_M,1);
 TN_old_passive_F.CorrByAnimal =  CorrelationsByAnimal(TN_old_passive_F,1);
 TN_young_passive2.CorrByAnimal = CorrelationsByAnimal(TN_young_passive2,1); 
 TN_old_passive.CorrByAnimal = CorrelationsByAnimal(TN_old_passive,1);
 %Compare2Anova(g1,g2,GroupNames,levels,SaveName,varnames,full_levels)

CorrDist.YoungAll = CorrsByDistance(TN_young_passive2,10,'YoungAll')
CorrDist.YoungM = CorrsByDistance(TN_young_passive_M,10,'YoungMale')
CorrDist.YoungF = CorrsByDistance(TN_young_passive_F,10,'YoungFemale')
CorrDist.AgingAll = CorrsByDistance(TN_old_passive,10,'AgingAll')
CorrDist.AgingM = CorrsByDistance(TN_old_passive_M,10,'AgingMale')
CorrDist.AgingF = CorrsByDistance(TN_old_passive_F,10,'AgingFemale')
 

CompareCorrsByDistance(CorrDist.YoungAll,CorrDist.AgingAll,'YoungAgingAll')
CompareCorrsByDistance(CorrDist.YoungM,CorrDist.AgingM,'YoungAgingMale')
CompareCorrsByDistance(CorrDist.YoungF,CorrDist.AgingF,'YoungAgingFemale')

 
 
 Corr.YoungMF.SignalCorr = CorrelationsAnalysisByCell(TN_young_passive_M.CorrByAnimal.LCorr',...
TN_young_passive_F.CorrByAnimal.LCorr',{'Male','Female'},{'Tones','20','10','0'},...
'SignalCorr-Cell') ;
                              
 Corr.YoungMF.NoiseCorr =  CorrelationsAnalysisByCell(TN_young_passive_M.CorrByAnimal.NCorr',...
TN_young_passive_F.CorrByAnimal.NCorr',{'Male','Female'},{'Tones','20','10','0'},...
'NoiseCorr-Cell') ;




 Corr.YoungAgingM.SignalCorr =  CorrelationsAnalysisByCell(TN_young_passive_M.CorrByAnimal.LCorr',...
TN_old_passive_M.CorrByAnimal.LCorr',{'Males-Young','Aging'},{'Tones','20','10','0'},...
'SignalCorr-Cell') ;

Corr.YoungAgingM.NoiseCorr =  CorrelationsAnalysisByCell(TN_young_passive_M.CorrByAnimal.NCorr',...
TN_old_passive_M.CorrByAnimal.NCorr',{'Males-Young','Aging'},{'Tones','20','10','0'},...
'NoiseCorr-Cell') ;

Corr.YoungAgingF.SignalCorr =  CorrelationsAnalysisByCell(TN_young_passive_F.CorrByAnimal.LCorr',...
TN_old_passive_F.CorrByAnimal.LCorr',{'Females-Young','Aging'},{'Tones','20','10','0'},...
'SignalCorr-Cell') ;

Corr.YoungAgingF.NoiseCorr =  CorrelationsAnalysisByCell(TN_young_passive_F.CorrByAnimal.NCorr',...
TN_old_passive_F.CorrByAnimal.NCorr',{'Females-Young','Aging'},{'Tones','20','10','0'},...
'NoiseCorr-Cell') ;


% Corr.AgingYoung.SignalCorr = CorrelationsAnalysisByCell(TN_young_passive2.CorrByAnimal.LCorr',...
% TN_old_passive2.CorrByAnimal.LCorr',{'Young','Aging'},{'Tones','20','10','0'},...
% 'SignalCorr-Cell') ;
% 
% Corr.AgingYoung.NoiseCorr =  CorrelationsAnalysisByCell(TN_young_passive2.CorrByAnimal.NCorr',...
% TN_old_passive2.CorrByAnimal.NCorr',{'Young','Aging'},{'Tones','20','10','0'},...
% 'NoiseCorr-Cell') ;



Corr.AgingYoung.SignalCorr = CorrelationsAnalysisByCell(TN_young_passive2.CorrByAnimal.LCorr',...
TN_old_passive.CorrByAnimal.LCorr',{'Young','Aging'},{'Tones','20','10','0'},...
'SignalCorr-Cell') ;

Corr.AgingYoung.NoiseCorr =  CorrelationsAnalysisByCell(TN_young_passive2.CorrByAnimal.NCorr',...
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
                                                                                                                                      
 
                           
  Corr.YoungAgingM.SignalCorrAnimal = CorrelationsAnalysisByAnimal(TN_young_passive_M.CorrByAnimal.LCorrAnimal',...
                                  TN_old_passive_M.CorrByAnimal.LCorrAnimal',...
                                  'SignalCorr-byAnimal',{'Males-Young','Aging'},...
                                  {'Tones','20','10','0'}) ;                           
                                                     
  
   Corr.YoungAgingM.NoiseCorrAnimal = CorrelationsAnalysisByAnimal(TN_young_passive_M.CorrByAnimal.NCorrAnimal',...
                                  TN_old_passive_M.CorrByAnimal.NCorrAnimal',...
                                  'NoiseCorr-byAnimal',{'Males-Young','Aging'},...
                                  {'Tones','20','10','0'}) ;
                              
                              
    Corr.YoungAgingF.SignalCorrAnimal = CorrelationsAnalysisByAnimal(TN_young_passive_F.CorrByAnimal.LCorrAnimal',...
                                  TN_old_passive_F.CorrByAnimal.LCorrAnimal',...
                                  'SignalCorr-byAnimal',{'Females-Young','Aging'},...
                                  {'Tones','20','10','0'}) ;                           
                                                     
  
   Corr.YoungAgingF.NoiseCorrAnimal = CorrelationsAnalysisByAnimal(TN_young_passive_F.CorrByAnimal.NCorrAnimal',...
                                  TN_old_passive_F.CorrByAnimal.NCorrAnimal',...
                                  'NoiseCorr-byAnimal',{'Females-Young','Aging'},...
                                  {'Tones','20','10','0'}) ;                             
                                                                                                                                           
                              
                              
  Corr.AgingYoung.SignalCorrAnimal = CorrelationsAnalysisByAnimal(TN_young_passive2.CorrByAnimal.LCorrAnimal',...
                                  TN_old_passive2.CorrByAnimal.LCorrAnimal',...
                                  'SignalCorr-byAnimal',{'Young','Aging'},...
                                  {'Tones','20','10','0'}) ;                           
  Corr.AgingYoung.NoiseCorrAnimal = CorrelationsAnalysisByAnimal(TN_young_passive2.CorrByAnimal.NCorrAnimal',...
                                  TN_old_passive2.CorrByAnimal.NCorrAnimal',...
                                  'NoiseCorr-byAnimal',{'Young','Aging'},...
                                  {'Tones','20','10','0'}) ; 
  
  Corr.AgingYoung.SignalCorrAnimal = CorrelationsAnalysisByAnimal(TN_young_passive2.CorrByAnimal.LCorrAnimal',...
                                  TN_old_passive.CorrByAnimal.LCorrAnimal',...
                                  'SignalCorr-byAnimal',{'Young','Aging'},...
                                  {'Tones','20','10','0'}) ;                           
  Corr.AgingYoung.NoiseCorrAnimal = CorrelationsAnalysisByAnimal(TN_young_passive2.CorrByAnimal.NCorrAnimal',...
                                  TN_old_passive.CorrByAnimal.NCorrAnimal',...
                                  'NoiseCorr-byAnimal',{'Young','Aging'},...
                                  {'Tones','20','10','0'}) ;   
 
                        



                     

     
                              

  % Plots that remove the effect of Group to focus on the effects of level
 Corr.AgingYoung.NoiseDiffCell = CorrelationsAnalysisByAnimal(TN_young_passive2.CorrByAnimal.NCorrDiff',...
                                  TN_old_passive.CorrByAnimal.NCorrDiff',...
                                  'NoiseDiff-byCell',{'Young','Aging-Rev'},...
                                  {'Tones','20','10','0'}) ;                           
                                                    
   Corr.AgingYoung.SignalDiffCell = CorrelationsAnalysisByAnimal(TN_young_passive2.CorrByAnimal.LCorrDiff',...
                                  TN_old_passive.CorrByAnimal.LCorrDiff',...
                                  'SignalDiff-byCell',{'Young','Aging-Rev'},...
                                  {'Tones','20','10','0'}) ;                           
                                                       

 Corr.AgingYoungM.NoiseDiffCell = CorrelationsAnalysisByAnimal(TN_young_passive_M.CorrByAnimal.NCorrDiff',...
                                  TN_old_passive_M.CorrByAnimal.NCorrDiff',...
                                  'NoiseDiff-byCell',{'Young','Aging-Rev'},...
                                  {'Tones','20','10','0'}) ;                           
                                                    
   Corr.AgingYoungM.SignalDiffCell = CorrelationsAnalysisByAnimal(TN_young_passive2.CorrByAnimal.LCorrDiff',...
                                  TN_old_passive_M.CorrByAnimal.LCorrDiff',...
                                  'SignalDiff-byCell',{'YoungM','AgingM'},...
                                  {'Tones','20','10','0'}) ;                           
                                                       





                              
  Corr.AgingYoung.NoiseDiffAnimal = CorrelationsAnalysisByAnimal(TN_young_passive2.CorrByAnimal.NCorrDiffAnimal',...
                                  TN_old_passive2.CorrByAnimal.NCorrDiffAnimal',...
                                  'NoiseDiff-byAnimal',{'Young','Aging-Rev'},...
                                  {'Tones','20','10','0'}) ;                           
                                                    
   Corr.AgingYoung.SignalDiffAnimal = CorrelationsAnalysisByAnimal(TN_young_passive2.CorrByAnimal.LCorrDiffAnimal',...
                                  TN_old_passive2.CorrByAnimal.LCorrDiffAnimal',...
                                  'SignalDiff-byAnimal',{'Young','Aging-Rev'},...
                                  {'Tones','20','10','0'}) ;                           
                                                    
 
             
  Corr.AgingYoungM.NoiseDiffAnimal = CorrelationsAnalysisByAnimal(TN_young_passive_M.CorrByAnimal.NCorrDiffAnimal',...
                                  TN_old_passive_M.CorrByAnimal.NCorrDiffAnimal',...
                                  'NoiseDiff-byAnimal',{'YoungM','AgingM'},...
                                  {'Tones','20','10','0'}) ;                           
                                                    
   Corr.AgingYoungM.SignalDiffAnimal = CorrelationsAnalysisByAnimal(TN_young_passive_M.CorrByAnimal.LCorrDiffAnimal',...
                                  TN_old_passive_M.CorrByAnimal.LCorrDiffAnimal',...
                                  'SignalDiff-byAnimal',{'YoungM','AgingM'},...
                                  {'Tones','20','10','0'}) ;                           
                              
 
    % Correlation by distance                    
                              
    Corr.                          
                              


%% Cluster comparision
All_Passive = Fluoro_to_Table([TN_young_passive_M.DataDirs,TN_young_passive_F.DataDirs,TN_old_passive_M.DataDirs,TN_old_passive_F.DataDirs] ,33000)
All_Passive.animalInfo = getAnimalInfo(All_Passive.DataDirs)
F_expts =  find(contains(All_Passive.animalInfo{:,2},'F'));
M_expts =  find(~contains(All_Passive.animalInfo{:,2},'F'));

old_expts = find(cellfun(@calyears,All_Passive.animalInfo{:,5}) > 0 );
young_expts = find(~ cellfun(@calyears,All_Passive.animalInfo{:,5}) > 0 );

F_idx = any(All_Passive.experiment_list == F_expts' ,2);
M_idx =  any(All_Passive.experiment_list == M_expts' ,2);
old_idx = any(All_Passive.experiment_list == old_expts' ,2);
young_idx = any(All_Passive.experiment_list == young_expts' ,2);


All_Passive.Clusters = Cluster_DF(All_Passive)
TN_old_passive.Combined_Classes =   All_Passive.Clusters.Clusters(old_idx);
TN_young_passive2.Combined_Classes = All_Passive.Clusters.Clusters(young_idx);

TN_old_passive_F.Combined_Classes =   All_Passive.Clusters.Clusters(old_idx & F_idx);
TN_old_passive_M.Combined_Classes =   All_Passive.Clusters.Clusters(old_idx & M_idx );
TN_young_passive_F.Combined_Classes = All_Passive.Clusters.Clusters( young_idx & F_idx);
TN_young_passive_M.Combined_Classes =  All_Passive.Clusters.Clusters(young_idx & M_idx );


Plot_Clusters(TN_young_passive2,TN_young_passive2.Combined_Classes,'Young_All')
Plot_Clusters(TN_old_passive,TN_old_passive.Combined_Classes,'Aging')

Plot_Clusters(TN_young_passive_M,TN_young_passive_M.Combined_Classes,'Young_Male')
Plot_Clusters(TN_young_passive_F,TN_young_passive_F.Combined_Classes,'Young_Female')

Plot_Clusters(TN_old_passive_M,TN_old_passive_M.Combined_Classes,'Aging_Male')
Plot_Clusters(TN_old_passive_F,TN_old_passive_F.Combined_Classes,'Aging_Female')



Diversity = struct()


Diversity.Stats.YoungAging = Diversityanalysis(TN_young_passive2,TN_old_passive,{'Young','Aging'});

Diversity.Stats.AgingYoungM =Diversityanalysis(TN_young_passive_M,TN_old_passive_M,{'Young','AgingM'});

Diversity.Stats.AgingYoungF = Diversityanalysis(TN_young_passive_F,TN_old_passive_F,{'YoungFemale','AgingFemale'});

Diversity.Stats.YoungMF = Diversityanalysis(TN_young_passive_M,TN_young_passive_F,{'YoungMale','Female'});



%  Bayes 
% TN_young_passive_M.Bayes = BayesClassiferPassive(TN_young_passive_M)
% TN_old_passive_M.Bayes = BayesClassiferPassive(TN_old_passive_M)
% TN_young_passive_F.Bayes = BayesClassiferPassive(TN_young_passive_F)
% TN_old_passive_F.Bayes = BayesClassiferPassive(TN_old_passive_F)

BayesStats.YoungAging = CompareBayes(TN_young_passive.Bayes,TN_old_passive_M.Bayes,...
    {'Young','Aging'},{'Tones','20','10','0'},...
'Bayes') ;


BayesStats.AgingYoungM = CompareBayes(TN_young_passive_M.Bayes,TN_old_passive_M.Bayes,...
      {'Males-Young','Aging'},{'Tones','20','10','0'},...
'Bayes') ;


BayesStats.YoungMF =  CompareBayes(TN_young_passive_M.Bayes,TN_young_passive_F.Bayes,...
        {'Males','Female'},{'Tones','20','10','0'},...
'Bayes') ;


BayesStats.AgingYoungF =  CompareBayes(TN_young_passive_F.Bayes,TN_old_passive_F.Bayes,...
        {'Females-Young','Aging'},{'Tones','20','10','0'},...
'Bayes') ;





