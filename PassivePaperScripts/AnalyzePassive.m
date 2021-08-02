

 %%  Intensity Comparision
 % Young MF
  [~,p1]=ttest2(PlotFluoroCDF(FRA_young_noise_F,'max',70),PlotFluoroCDF(FRA_young_noise_M,'max',70)) ; 
  [~,p2]=ttest2(PlotFluoroCDF(FRA_young_noise_F,'mean'),PlotFluoroCDF(FRA_young_noise_M,'mean')) 

% Young V Old
 [~,p3]=ttest2(PlotFluoroCDF(FRA_young_noise2,'max'),PlotFluoroCDF(FRA_old_noise,'max')) 
  [~,p4]=ttest2(PlotFluoroCDF(FRA_young_noise2,'mean'),PlotFluoroCDF(FRA_old_noise,'mean')) 
  
  
  IntesnityStats = struct()
  IntensityStats.MaleFemale_Young.Max = p1;
  IntensityStats.MaleFemale_Young.Mean = p2;
  IntensityStats.OldYoung.Max = p3;
  IntensityStatsOldYoung.Mean = p4;

clear p1 p2 p3 p4
  
 %% Temporal Analysis
     % Timing Noise
        % young MF
[~,young_timing_noise_F]=TemporalAnalysisByAnimal(FRA_young_noise_F,70,'Timing_Young_F_noise')
    
[~,young_timing_noise_M]=TemporalAnalysisByAnimal(FRA_young_noise_M,70,'Timing_Young_M_noise')

 mfn =CompareOldYoung(young_timing_noise_F',young_timing_noise_M',[1:size(young_timing_noise_F ,2)],'CompareTiming_YoungFemaleMale_Noise')


         % Old vs young
[~,old_timing_noise]=TemporalAnalysisByAnimal(FRA_old_noise,70,'Timing_Old_noise');
    
[~,young_timing_noise]=TemporalAnalysisByAnimal(FRA_young_noise2,70,'Timing_Young_noise');

oyn=CompareOldYoung(old_timing_noise',young_timing_noise',[1:size(young_timing_noise,2)],'CompareTiming_OldYoungTiming-Noise')


    % Timing_quiet 
[~,young_timing_quiet_F]=TemporalAnalysisByAnimal(FRA_young_noise_F,inf,'Timing_Young_F_quiet');
    
[~,young_timing_quiet_M]=TemporalAnalysisByAnimal(FRA_young_noise_M,inf,'Timing_Young_M_quiet');

mfq=CompareOldYoung(young_timing_quiet_F',young_timing_quiet_M',[1:size(young_timing_noise_M,2)],'CompareTiming_YoungFemaleMale_Quiet')



[~,old_timing_quiet]=TemporalAnalysisByAnimal(FRA_old_noise,inf,'Timing_Old_quiet');
    
[~,young_timing_quiet]=TemporalAnalysisByAnimal(FRA_young_noise2,inf,'Timing_Young_quiet');

oyq=CompareOldYoung(old_timing_quiet',young_timing_quiet',[1:size(young_timing_quiet,2)],'CompareTiming_OldYoung_Quiet')

TimingStats = struct()

TimingStats.OldYoung.Noise = oyn;
TimingStats.OldYoung.Quiet = oyq;

TimingStats.Young_FemaleMale.Noise = mfn;
TimingStats.Young_FemaleMale.Quiet = mfq;

save('TimingStats.mat','TimingStats')

clear oyn oyq mfq mfn

clear young_timing_noise young_timing_quiet 
clear old_timing_noise old_timing_quiet
clear young_timing_noise_M young_timing_noise_F
clear young_timing_quiet_M young_timing_quiet_F




%% Bandwidth Comparision
BandwithAll = struct()

BandwidthAll.Young.All.Pos = BandwidthAnalysis(FRA_young_noise2, 'BRFS',0,'Pos',.5,'YoungAll-Pos')
BandwidthAll.Young.All.Neg =BandwidthAnalysis(FRA_young_noise2, 'BRFS',0,'Neg',.5,'YoungAll-Neg')

BandwidthAll.Old.All.Pos = BandwidthAnalysis(FRA_old_noise, 'BRFS',0,'Pos',.5,'OldAll-Pos')
BandwidthAll.Old.All.Neg =BandwidthAnalysis(FRA_old_noise, 'BRFS',0,'Neg',.5,'OldAll-Neg')

BandwidthAll.Young.F.Pos = BandwidthAnalysis(FRA_young_noise_F, 'BRFS',0,'Pos',.5,'YoungF-Pos')
BandwidthAll.Young.F.Neg =BandwidthAnalysis(FRA_young_noise_F, 'BRFS',0,'Neg',.5,'YoungF-Neg')

BandwidthAll.Young.M.Pos = BandwidthAnalysis(FRA_young_noise_M, 'BRFS',0,'Pos',.5,'YoungM-Pos')
BandwidthAll.Young.M.Neg =BandwidthAnalysis(FRA_young_noise_M, 'BRFS',0,'Neg',.5,'YoungM-Neg')



% Bandwidth Stats = 
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

                             
      % StatsByAnimal                       
BandwithAll.Stats.OldYoung.NegByAnimal = CompareOldYoung(BandwidthAll.Old.All.Neg{4},...
                                 BandwidthAll.Young.All.Neg{4},...
                                  {'Tones','20','10','0'},...
                                 'OldYoungAll-NegByAnimal');
                             
BandwithAll.Stats.OldYoung.PosByAnimal = CompareOldYoung(BandwidthAll.Old.All.Pos{4}',...
                                 BandwidthAll.Young.All.Pos{4}',...
                                  {'Tones','20','10','0'},...
                                 'OldYoungAll-PosByAnimal');                             
                             
                             
Bandwidth.Stats.YoungMF.PosByAnimal = CompareOldYoung(BandwidthAll.Young.F.Pos{4},...
                                 BandwidthAll.Young.M.Pos{4},...
                                  {'Tones','20','10','0'},...
                                 'YoungMF-PosByAnimal');   
                             
Bandwidth.Stats.YoungMF.NegByAnimal = CompareOldYoung(BandwidthAll.Young.F.Neg{4},...
                                 BandwidthAll.Young.M.Neg{4},...
                                  {'Tones','20','10','0'},...
                                 'YoungMF-NegByAnimal'); 

                                                        
                            



%
