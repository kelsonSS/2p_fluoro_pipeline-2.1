% PassiveCompareSexScript

% generate datasets seperated by sex
 old_F_idx = contains(FRA_old_noise2.AnimalInfo{:,2},'F'); 
young_F_idx= contains(FRA_young_noise2.AnimalInfo{:,2},'F');

 
 FRA_old_noise_F = Fluoro_to_Table(FRA_old_noise2.DataDirs(old_F_idx),33000);
  FRA_old_noise_M = Fluoro_to_Table(FRA_old_noise2.DataDirs(~old_F_idx),33000);

  
   FRA_young_noise_F = Fluoro_to_Table(FRA_young_noise2.DataDirs(young_F_idx),33000);
  FRA_young_noise_M = Fluoro_to_Table(FRA_young_noise2.DataDirs(~young_F_idx),33000);
  
  % Age Intensity
  [~,p]=ttest2(PlotFluoroCDF(FRA_old_noise,'mean'),PlotFluoroCDF(FRA_young_noise2,'mean')) 
  [~,p]=ttest2(PlotFluoroCDF(FRA_old_noise,'max'),PlotFluoroCDF(FRA_young_noise2,'max')) 
  % young Intensity
  [~,p]=ttest2(PlotFluoroCDF(FRA_young_noise_F,'max'),PlotFluoroCDF(FRA_young_noise_M,'max')) 
  [~,p]=ttest2(PlotFluoroCDF(FRA_young_noise_F,'mean'),PlotFluoroCDF(FRA_young_noise_M,'mean')) 
  % old Intensity
  
    [~,p]=ttest2(PlotFluoroCDF(FRA_old_noise_F,'mean'),PlotFluoroCDF(FRA_old_noise_M,'mean'))
     [~,p]=ttest2(PlotFluoroCDF(FRA_old_noise_F,'max'),PlotFluoroCDF(FRA_old_noise_M,'max'))
    
    % Timing Noise
[~,young_timing_noise_F]=TemporalAnalysisByAnimal(FRA_young_noise_F,70)
    
[~,young_timing_noise_M]=TemporalAnalysisByAnimal(FRA_young_noise_M,70)

CompareOldYoung(young_timing_noise_F',young_timing_noise_M',[1:size(young_timing_noise_M,2)])



[~,old_timing_noise_F]=TemporalAnalysisByAnimal(FRA_old_noise_F,70);
    
[~,old_timing_noise_M]=TemporalAnalysisByAnimal(FRA_old_noise_M,70);

CompareOldYoung(old_timing_noise_F',old_timing_noise_M',[1:size(old_timing_noise_M,2)])


% Timing_quiet 
[~,young_timing_noise_F]=TemporalAnalysisByAnimal(FRA_young_noise_F,inf)
    
[~,young_timing_noise_M]=TemporalAnalysisByAnimal(FRA_young_noise_M,inf)

CompareOldYoung(young_timing_noise_F',young_timing_noise_M',[1:size(young_timing_noise_M,2)])



[~,old_timing_noise_F]=TemporalAnalysisByAnimal(FRA_old_noise_F,inf);
    
[~,old_timing_noise_M]=TemporalAnalysisByAnimal(FRA_old_noise_M,inf);

CompareOldYoung(old_timing_noise_F',old_timing_noise_M',[1:size(old_timing_noise_M,2)])



%Bandwidth Comparision
BandwithAll = struct()

BandwidthAll.Young.All.Pos = BandwidthAnalysis(FRA_young_noise2, 'BRFS',0,'Pos',.5,'YoungAll-Pos')
BandwidthAll.Young.All.Neg =BandwidthAnalysis(FRA_young_noise2, 'BRFS',0,'Neg',.5,'YoungAll-Neg')

BandwidthAll.Old.All.Pos = BandwidthAnalysis(FRA_old_noise2, 'BRFS',0,'Pos',.5,'OldAll-Pos')
BandwidthAll.Old.All.Neg =BandwidthAnalysis(FRA_old_noise2, 'BRFS',0,'Neg',.5,'OldAll-Neg')

BandwidthAll.Young.F.Pos = BandwidthAnalysis(FRA_young_noise_F, 'BRFS',0,'Pos',.5,'YoungF-Pos')
BandwidthAll.Young.F.Neg =BandwidthAnalysis(FRA_young_noise_F, 'BRFS',0,'Neg',.5,'YoungF-Neg')

BandwidthAll.Young.M.Pos = BandwidthAnalysis(FRA_young_noise_M, 'BRFS',0,'Pos',.5,'YoungM-Pos')
BandwidthAll.Young.M.Neg =BandwidthAnalysis(FRA_young_noise_M, 'BRFS',0,'Neg',.5,'YoungM-Neg')


BandwidthAll.Old.F.Pos = BandwidthAnalysis(FRA_old_noise_F, 'BRFS',0,'Pos',.5,'OldF-Pos')
BandwidthAll.Old.F.Neg =BandwidthAnalysis(FRA_old_noise_F, 'BRFS',0,'Neg',.5,'OldF-Neg')

BandwidthAll.Old.M.Pos = BandwidthAnalysis(FRA_old_noise_M, 'BRFS',0,'Pos',.5,'OldM-Pos')
BandwidthAll.Old.M.Neg =BandwidthAnalysis(FRA_old_noise_M, 'BRFS',0,'Neg',.5,'OldM-Neg')


% Stats = 
BandwithAll.Stats.OldYoung.Neg = CompareOldYoung(BandwidthAll.Old.All.Neg{1},...
                                 BandwidthAll.Young.All.Neg{1},...
                                  {'Tones','20','10','0'},...
                                 'OldYoungAll-Neg');
                             
BandwithAll.Stats.OldYoung.Pos = CompareOldYoung(BandwidthAll.Old.All.Pos{1},...
                                 BandwidthAll.Young.All.Pos{1},...
                                  {'Tones','20','10','0'},...
                                 'OldYoungAll-Pos');                             
                             
                             
Bandwidth.Stats.YoungMF.Pos = CompareOldYoung(BandwidthAll.Young.F.Pos{1},...
                                 BandwidthAll.Young.M.Pos{1},...
                                  {'Tones','20','10','0'},...
                                 'YoungMF-Pos');   
                             
Bandwidth.Stats.YoungMF.Neg = CompareOldYoung(BandwidthAll.Young.F.Neg{1},...
                                 BandwidthAll.Young.M.Neg{1},...
                                  {'Tones','20','10','0'},...
                                 'YoungMF-Neg');                                
                             




% Correlations 

 FRA_young_noise_F.CorrByAnimal = CorrelationsByAnimal(FRA_young_noise_F,1);
 FRA_young_noise_M.CorrByAnimal = CorrelationsByAnimal(FRA_young_noise_M,1);
 
 Corr.YoungMF.SignalCorr = CompareOldYoung(FRA_young_noise_F.CorrByAnimal.LCorrAnimal,...
                                  FRA_young_noise_M.CorrByAnimal.LCorrAnimal,...
                                  {'Tones','20','10','0'}) ;
                              
                              
 Corr.YoungMF.NoiseCorr = CompareOldYoung(FRA_young_noise_F.CorrByAnimal.NCorrAnimal,...
                                  FRA_young_noise_M.CorrByAnimal.NCorrAnimal,...
                                  {'Tones','20','10','0'}) ;
                                                                        
                                          


% Cluster comparision
ALL_Passive = Fluoro_to_Table([FRA_young_noise2.DataDirs,FRA_old_noise.DataDirs] ,33000)
F_expts =  find(contains(All_Passive.animalInfo{:,2},'F'));
old_expts = find(cellfun(@calyears,All_Passive.animalInfo{:,4}) > 0 );

F_idx = any(All_Passive.experiment_list == F_expts' ,2);
old_idx = any(All_Passive.experiment_list == old_expts' ,2);


All_Passive.Clusters = Cluster_DF(All_Passive)
FRA_old_noise2.Combined_Classes =   All_Passive.Clusters.Clusters(old_idx);
FRA_young_noise2.Combined_Classes = All_Passive.Clusters.Clusters(~old_idx);

FRA_old_noise_F.Combined_Classes =   All_Passive.Clusters.Clusters(old_idx & F_idx);
FRA_old_noise_M.Combined_Classes =   All_Passive.Clusters.Clusters(old_idx & ~F_idx);
FRA_young_noise_F.Combined_Classes = All_Passive.Clusters.Clusters(~old_idx & F_idx);
FRA_young_noise_M.Combined_Classes =  All_Passive.Clusters.Clusters(~old_idx & ~F_idx);

Plot_ClusterDiversity(FRA_old_noise2.Combined_Classes)
Plot_ClusterDiversity(FRA_o_noise_M.Combined_Classes)
Plot_Clusters(FRA_young_noise,FRA_young_noise.Combined_Classes)




