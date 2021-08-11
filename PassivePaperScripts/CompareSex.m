% PassiveCompareSexScript

% get all old animals less than 20 months old.
FRA_old_noise2 =Fluoro_to_Table( FRA_old_noise_All.DataDirs( ...
      cellfun(@calmonths,FRA_old_noise_All.AnimalInfo.Months_old) < 20),33000)
  FRA_old_noise2.AnimalInfo = getAnimalInfo(FRA_old_noise2.DataDirs)
% generate datasets seperated by sex
 old_F_idx = contains(FRA_old_noise2.AnimalInfo{:,2},'F'); 
young_F_idx= contains(FRA_young_noise2.AnimalInfo{:,2},'F');

 
 FRA_old_noise_F = Fluoro_to_Table(FRA_old_noise2.DataDirs(old_F_idx),33000);
  FRA_old_noise_M = Fluoro_to_Table(FRA_old_noise2.DataDirs(~old_F_idx),33000);

  
   FRA_young_noise2_F = Fluoro_to_Table(FRA_young_noise2.DataDirs(young_F_idx),33000);
  FRA_young_noise2_M = Fluoro_to_Table(FRA_young_noise2.DataDirs(~young_F_idx),33000);
  

  
  %% Intensity
  
  % Age Intensity
  Intensity = struct()
  
  Intensity.YoungAging.Quiet.Mean =  FluoroIntensityAnalysis(FRA_young_noise2,FRA_old_noise2,...
      'mean',{'Young','Aging'},Inf)
  
   Intensity.YoungAging.Noise.Mean =  FluoroIntensityAnalysis(FRA_young_noise2,FRA_old_noise2,...
      'mean',{'Young','Aging'},70)
  
  Intensity.YoungAging.Quiet.Max =  FluoroIntensityAnalysis(FRA_young_noise2,FRA_old_noise2,...
      'max',{'Young','Aging'},Inf)
  
   Intensity.YoungAging.Noise.Max =  FluoroIntensityAnalysis(FRA_young_noise2,FRA_old_noise2,...
      'max',{'Young','Aging'},70)
  
  
  % Sex 
  
    Intensity.YoungMF.Quiet.Mean = FluoroIntensityAnalysis(FRA_young_noise_M,FRA_young_noise_F,...
      'mean',{'Young-M','Young-F'},Inf)
  
   Intensity.YoungMF.Quiet.Max =  FluoroIntensityAnalysis(FRA_young_noise_M,FRA_young_noise_F,...
      'max',{'Young-M','Young-F'},Inf)
  
    Intensity.YoungMF.Noise.Mean = FluoroIntensityAnalysis(FRA_young_noise_M,FRA_young_noise_F,...
      'mean',{'Young-M','Young-F'},70)
  
   Intensity.YoungMF.Noise.Max =  FluoroIntensityAnalysis(FRA_young_noise_M,FRA_young_noise_F,...
      'max',{'Young-M','Young-F'},70)
   
  % old M vs young M
  
  
   Intensity.OldYoungM.Quiet.Mean = FluoroIntensityAnalysis(FRA_old_noise_M,FRA_young_noise_M,...
      'mean',{'Old-M','Young-M'},Inf)
  
   Intensity.OldYoungM.Quiet.Max =  FluoroIntensityAnalysis(FRA_old_noise_M,FRA_young_noise_M,...
      'max',{'Old-M','Young-M'},Inf)
  
    Intensity.OldYoungM.Noise.Mean = FluoroIntensityAnalysis(FRA_old_noise_M,FRA_young_noise_M,...
      'mean',{'Old-M','Young-M'},70)
  
   Intensity.OldYoungM.Noise.Max =  FluoroIntensityAnalysis(FRA_old_noise_M,FRA_young_noise_M,...
      'max',{'Old-M','Young-M'},70)
  
    % old F vs young F
  
  
   Intensity.OldYoungF.Quiet.Mean = FluoroIntensityAnalysis(FRA_young_noise_F,FRA_old_noise_F,...
      'mean',{'Young-F','Old-F'},Inf)
  
   Intensity.OldYoungF.Quiet.Max =  FluoroIntensityAnalysis(FRA_young_noise_F,FRA_old_noise_F,...
      'max',{'Young-F','Old-F'},Inf)
  
    Intensity.OldYoungF.Noise.Mean = FluoroIntensityAnalysis(FRA_young_noise_F,FRA_old_noise_F,...
      'mean',{'Young-F','Old-F'},70)
  
   Intensity.OldYoungF.Noise.Max =  FluoroIntensityAnalysis(FRA_young_noise_F,FRA_old_noise_F,...
      'max',{'Young-F','Old-F'},70)
  
  
  
save 'IntensityStats.mat' 'Intensity' 

  

  %% Timing
    %% Timing Noise 
      % Young
[~,young_timing_noise_F]=TemporalAnalysisByAnimal(FRA_young_noise_F,70)
    
[~,young_timing_noise_M]=TemporalAnalysisByAnimal(FRA_young_noise_M,70)

[~,old_timing_noise_F]=TemporalAnalysisByAnimal(FRA_old_noise_F,70,'Timing_Old_F_noise');

[~,old_timing_noise_M]=TemporalAnalysisByAnimal(FRA_old_noise_M,70);

[~,old_timing_noise]=TemporalAnalysisByAnimal(FRA_old_noise2,70);
    
[~,young_timing_noise]=TemporalAnalysisByAnimal(FRA_young_noise2,70);


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
      

%% Timing_quiet 

% 
[~,young_timing_noise_F]=TemporalAnalysisByAnimal(FRA_young_noise_F,inf)
    
[~,young_timing_noise_M]=TemporalAnalysisByAnimal(FRA_young_noise_M,inf)

CompareOldYoung(young_timing_noise_F',young_timing_noise_M',[1:size(young_timing_noise_M,2)])




[~,old_timing_noise_F]=TemporalAnalysisByAnimal(FRA_old_noise_F,inf,'Timing_Old_F_quiet');
    
[~,old_timing_noise_M]=TemporalAnalysisByAnimal(FRA_old_noise_M,inf);

CompareOldYoung(old_timing_noise_F',old_timing_noise_M',[1:size(old_timing_noise_M,2)])

% Old_M vs Young_M
      CompareOldYoung(young_timing_noise_M',old_timing_noise_M',[1:size(old_timing_noise_M,2)])







%% Bandwidth Comparision
BandwithAll = struct()

% Using BRFS 
BandwidthAll.Young.All.Pos = BandwidthAnalysis(FRA_young_noise2, 'BRFS',0,'Pos',.5,'YoungAll-Pos')
BandwidthAll.Young.All.Neg = BandwidthAnalysis(FRA_young_noise2, 'BRFS',0,'Neg',.5,'YoungAll-Neg')

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


% using traditional_BW
BandwidthAll.Young.All.BW_Pos = BandwidthAnalysis(FRA_young_noise2, 'interp',0,'Pos',.5,'YoungAll-BW-Pos')
BandwidthAll.Young.All.BW_Neg =BandwidthAnalysis(FRA_young_noise2, 'interp',0,'Neg',.5,'YoungAll-BW-Neg')

BandwidthAll.Old.All.BW_Pos = BandwidthAnalysis(FRA_old_noise2, 'interp',0,'Pos',.5,'OldAll-BW-Pos')
BandwidthAll.Old.All.BW_Neg =BandwidthAnalysis(FRA_old_noise2, 'interp',0,'Neg',.5,'OldAll-BW-Neg')

BandwidthAll.Young.F.BW_Pos = BandwidthAnalysis(FRA_young_noise_F, 'interp',0,'Pos',.5,'YoungF-BW_Pos')
BandwidthAll.Young.F.BW_Neg =BandwidthAnalysis(FRA_young_noise_F, 'interp',0,'Neg',.5,'YoungF-BW_Neg')

BandwidthAll.Young.M.BW_Pos = BandwidthAnalysis(FRA_young_noise_M, 'interp',0,'Pos',.5,'YoungM-BW_Pos')
BandwidthAll.Young.M.BW_Neg =BandwidthAnalysis(FRA_young_noise_M, 'interp',0,'Neg',.5,'YoungM-BW_Neg')

BandwidthAll.Old.M.BW_Pos = BandwidthAnalysis(FRA_old_noise_M, 'interp',0,'Pos',.5,'OldM-BW_Pos')
BandwidthAll.Old.M.BW_Neg = BandwidthAnalysis(FRA_old_noise_M, 'interp',0,'Neg',.5,'OldM-BW_Neg')


BandwidthAll.Old.F.BW_Pos = BandwidthAnalysis(FRA_old_noise_F, 'interp',0,'Pos',.5,'OldF-BW_Pos')
BandwidthAll.Old.F.BW_Neg =BandwidthAnalysis(FRA_old_noise_F, 'interp',0,'Neg',.5,'OldF-BW_Neg')



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


BandwidthAll.Stats.OldYoung.Neg = CompareOldYoung(BandwidthAll.Old.All.Neg,...
                                 BandwidthAll.Young.All.Neg,{'Aging','Young'},...
                                  {'Tones','20','10','0'},...
                                 'OldYoungAll-Neg');
                             
BandwidthAll.Stats.OldYoung.Pos = CompareOldYoung(BandwidthAll.Old.All.Pos,...
                                 BandwidthAll.Young.All.Pos,{'Aging','Young'},......
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
               
BandwidthAll.Stats.YoungOldM.Pos = CompareOldYoung(BandwidthAll.Old.M.Pos,...
                                 BandwidthAll.Young.M.Pos,{'Old-M','Young-M'},...
                                  {'Tones','20','10','0'},...
                                 'YoungOldM-Pos');                               
                             
BandwidthAll.Stats.YoungOldM.Neg = CompareOldYoung(BandwidthAll.Old.M.Neg,...
                                 BandwidthAll.Young.M.Neg,{'Old-M','Young-M'},...
                                  {'Tones','20','10','0'},...
                                 'YoungOldM-Neg');     
                                                          
     % Traditional BW

BandwidthAll.Stats.OldYoung.BW_Pos = CompareOldYoung(BandwidthAll.Old.All.BW_Pos,...
                                 BandwidthAll.Young.All.BW_Pos,{'Aging','Young'},......
                                  {'Tones','20','10','0'},...
                                 'OldYoungAll-BW_Pos');
            

BandwidthAll.Stats.OldYoung.BW_Neg = CompareOldYoung(BandwidthAll.Old.All.BW_Neg,...
                                 BandwidthAll.Young.All.BW_Neg,{'Aging','Young'},......
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
                              
                     
                                                               
  BandwidthAll.Stats.YoungOldF.BW_Pos = CompareOldYoung(BandwidthAll.Old.F.BW_Pos,...
                                 BandwidthAll.Young.F.BW_Pos,{'Old-F','Young-F'},...
                                  {'Tones','20','10','0'},...
                                 'YoungOldF-BW_Pos');  
                             
                              
 BandwidthAll.Stats.YoungOldF.BW_Neg = CompareOldYoung(BandwidthAll.Old.F.BW_Neg,...
                                 BandwidthAll.Young.F.BW_Neg,{'Old-F','Young-F'},...
                                  {'Tones','20','10','0'},...
                                 'YoungOldF-BW_Neg');  
       
                                     
                            
                             
                             




%% Correlations 

 FRA_young_noise_F.CorrByAnimal = CorrelationsByAnimal(FRA_young_noise_F,1);
 FRA_young_noise_M.CorrByAnimal = CorrelationsByAnimal(FRA_young_noise_M,1);
 FRA_old_noise_M.CorrByAnimal = CorrelationsByAnimal(FRA_old_noise_M,1);
 
 %Compare2Anova(g1,g2,GroupNames,levels,SaveName,varnames,full_levels)

 
 
 
 Corr.YoungMF.SignalCorr = Compare2Anova(FRA_young_noise_M.CorrByAnimal.LCorr',...
FRA_young_noise_F.CorrByAnimal.LCorr',{'Male','Female'},{'Tones','20','10','0'},...
'SignalCorr-Cell') ;
                              
 Corr.YoungMF.NoiseCorr =  Compare2Anova(FRA_young_noise_M.CorrByAnimal.NCorr',...
FRA_young_noise_F.CorrByAnimal.NCorr',{'Male','Female'},{'Tones','20','10','0'},...
'NoiseCorr-Cell') ;



 Corr.YoungOldM.SignalCorr =  Compare2Anova(FRA_young_noise_M.CorrByAnimal.LCorr',...
FRA_old_noise_M.CorrByAnimal.LCorr',{'Males-Young','Old'},{'Tones','20','10','0'},...
'SignalCorr-Cell') ;

Corr.YoungOldM.NoiseCorr =  Compare2Anova(FRA_young_noise_M.CorrByAnimal.NCorr',...
FRA_old_noise_M.CorrByAnimal.NCorr',{'Males-Young','Old'},{'Tones','20','10','0'},...
'NoiseCorr-Cell') ;

Corr.YoungOldF.SignalCorr =  Compare2Anova(FRA_young_noise_F.CorrByAnimal.LCorr',...
FRA_old_noise_F.CorrByAnimal.LCorr',{'Females-Young','Old'},{'Tones','20','10','0'},...
'SignalCorr-Cell') ;

Corr.YoungOldF.NoiseCorr =  Compare2Anova(FRA_young_noise_F.CorrByAnimal.NCorr',...
FRA_old_noise_F.CorrByAnimal.NCorr',{'Females-Young','Old'},{'Tones','20','10','0'},...
'NoiseCorr-Cell') ;


Corr.OldYoung.SignalCorr = Compare2Anova(FRA_young_noise2.CorrByAnimal.LCorr',...
FRA_old_noise2.CorrByAnimal.LCorr',{'Young','Aging'},{'Tones','20','10','0'},...
'SignalCorr-Cell') ;

Corr.OldYoung.NoiseCorr =  Compare2Anova(FRA_young_noise2.CorrByAnimal.NCorr',...
FRA_old_noise2.CorrByAnimal.NCorr',{'Young','Aging'},{'Tones','20','10','0'},...
'NoiseCorr-Cell') ;


 
  Corr.YoungMF.SignalCorrAnimal = CorrelationsAnalysisByAnimal(FRA_young_noise_M.CorrByAnimal.LCorrAnimal',...
                                  FRA_young_noise_F.CorrByAnimal.LCorrAnimal',...
                                  'SignalCorr-byAnimal',{'Male','Female'},...
                                  {'Tones','20','10','0'}) ;  
                              
 Corr.YoungMF.NoiseCorrAnimal = CorrelationsAnalysisByAnimal(FRA_young_noise_M.CorrByAnimal.NCorrAnimal',...
                                  FRA_young_noise_F.CorrByAnimal.NCorrAnimal',...
                                  'NoiseCorr-byAnimal',{'Male','Female'},...
                                  {'Tones','20','10','0'}) ;                           
                                                                                                                                      
 
                           
  Corr.YoungOldM.SignalCorrAnimal = CorrelationsAnalysisByAnimal(FRA_young_noise_M.CorrByAnimal.LCorrAnimal',...
                                  FRA_old_noise_M.CorrByAnimal.LCorrAnimal',...
                                  'SignalCorr-byAnimal',{'Males-Young','Aging'},...
                                  {'Tones','20','10','0'}) ;                           
                                                     
  
   Corr.YoungOldM.NoiseCorrAnimal = CorrelationsAnalysisByAnimal(FRA_young_noise_M.CorrByAnimal.NCorrAnimal',...
                                  FRA_old_noise_M.CorrByAnimal.NCorrAnimal',...
                                  'NoiseCorr-byAnimal',{'Males-Young','Aging'},...
                                  {'Tones','20','10','0'}) ;
                              
                              
    Corr.OldYoungF.SignalCorrAnimal = CorrelationsAnalysisByAnimal(FRA_young_noise_F.CorrByAnimal.LCorrAnimal',...
                                  FRA_old_noise_F.CorrByAnimal.LCorrAnimal',...
                                  'SignalCorr-byAnimal',{'Females-Young','Aging'},...
                                  {'Tones','20','10','0'}) ;                           
                                                     
  
   Corr.OldYoungF.NoiseCorrAnimal = CorrelationsAnalysisByAnimal(FRA_young_noise_F.CorrByAnimal.NCorrAnimal',...
                                  FRA_old_noise_F.CorrByAnimal.NCorrAnimal',...
                                  'NoiseCorr-byAnimal',{'Females-Young','Aging'},...
                                  {'Tones','20','10','0'}) ;                             
                                                                                                                                           
                              
                              
  Corr.OldYoung.SignalCorrAnimal = CorrelationsAnalysisByAnimal(FRA_young_noise2.CorrByAnimal.LCorrAnimal',...
                                  FRA_old_noise2.CorrByAnimal.LCorrAnimal',...
                                  'SignalCorr-byAnimal',{'Young','Aging'},...
                                  {'Tones','20','10','0'}) ;                           
  Corr.OldYoung.NoiseCorrAnimal = CorrelationsAnalysisByAnimal(FRA_young_noise2.CorrByAnimal.NCorrAnimal',...
                                  FRA_old_noise2.CorrByAnimal.NCorrAnimal',...
                                  'NoiselCorr-byAnimal',{'Young','Aging'},...
                                  {'Tones','20','10','0'}) ;   
                        



                     

     
                              

  % Plots that remove the effect of Group to focus on the effects of level
   

                              
  Corr.OldYoung.NoiseDiffAnimal = CorrelationsAnalysisByAnimal(FRA_young_noise2.CorrByAnimal.NCorrDiffAnimal',...
                                  FRA_old_noise2.CorrByAnimal.NCorrDiffAnimal',...
                                  'NoiseDiff-byAnimal',{'Young','Aging'},...
                                  {'Tones','20','10','0'}) ;                           
                                                    
   Corr.OldYoung.SignalDiffAnimal = CorrelationsAnalysisByAnimal(FRA_young_noise2.CorrByAnimal.LCorrDiffAnimal',...
                                  FRA_old_noise2.CorrByAnimal.LCorrDiffAnimal',...
                                  'SignalDiff-byAnimal',{'Young','Aging'},...
                                  {'Tones','20','10','0'}) ;                           
                                                    
 
                              
                              
 
                       
                              
                              
                              


%% Cluster comparision
All_Passive = Fluoro_to_Table([FRA_young_noise_M.DataDirs,FRA_young_noise_F.DataDirs,FRA_old_noise_M.DataDirs,FRA_old_noise_F.DataDirs] ,33000)
All_Passive.animalInfo = getAnimalInfo(All_Passive.DataDirs)
F_expts =  find(contains(All_Passive.animalInfo{:,2},'F'));
M_expts =  find(~contains(All_Passive.animalInfo{:,2},'F'));

old_expts = find(cellfun(@calyears,All_Passive.animalInfo{:,5}) > 0 );

F_idx = any(All_Passive.experiment_list == F_expts' ,2);
M_idx =  any(All_Passive.experiment_list == M_expts' ,2);
old_idx = any(All_Passive.experiment_list == old_expts' ,2);


All_Passive.Clusters = Cluster_DF(All_Passive)
FRA_old_noise2.Combined_Classes =   All_Passive.Clusters.Clusters(old_idx);
FRA_young_noise2.Combined_Classes = All_Passive.Clusters.Clusters(~old_idx);

FRA_old_noise_F.Combined_Classes =   All_Passive.Clusters.Clusters(old_idx & F_idx);
FRA_old_noise_M.Combined_Classes =   All_Passive.Clusters.Clusters(old_idx & M_idx );
FRA_young_noise_F.Combined_Classes = All_Passive.Clusters.Clusters((~old_idx) & F_idx);
FRA_young_noise_M.Combined_Classes =  All_Passive.Clusters.Clusters((~old_idx) & M_idx );

Plot_ClusterDiversity(FRA_old_noise2.Combined_Classes)
Plot_ClusterDiversity(FRA_young_noise_M.Combined_Classes)
title('Young M Cluster Diversity')
Plot_ClusterDiversity(FRA_young_noise_F.Combined_Classes)
title('Young F Cluster Diversity')
Plot_ClusterDiversity(FRA_old_noise_M.Combined_Classes)
title('Old M Cluster Diversity')
Plot_Clusters(FRA_young_noise2,FRA_young_noise2.Combined_Classes)
title('Old F Cluster Diversity')
Plot_Clusters(FRA_old_noise_F,FRA_old_noise_F.Combined_Classes)



Diversity = struct()


Diversity.Stats.YoungOld = Diversityanalysis(FRA_young_noise2,FRA_old_noise2,{'Young','Aging'})

Diversity.Stats.OldYoungM =Diversityanalysis(FRA_young_noise_M,FRA_old_noise_M,{'Young','OldM'});

Diversity.Stats.YoungMF = Diversityanalysis(FRA_young_noise_M,FRA_young_noise_F,{'YoungMale','Female'});

Diversity.Stats.OldYoungF = Diversityanalysis(FRA_young_noise_F,FRA_old_noise_F,{'YoungFemale','OldFemale'});


%  Bayes 
FRA_young_noise_M.Bayes = BayesClassiferPassive(FRA_young_noise_M)
FRA_old_noise_M.Bayes = BayesClassiferPassive(FRA_old_noise_M)
FRA_young_noise_F.Bayes = BayesClassiferPassive(FRA_young_noise_F)
FRA_old_noise_F.Bayes = BayesClassiferPassive(FRA_old_noise_F)

BayesStats.YoungOld = CompareBayes(FRA_young_noise.Bayes,FRA_old_noise_M.Bayes,...
    {'Young','Aging'},{'Tones','20','10','0'},...
'Bayes') ;

BayesStats.OldYoungM = CompareBayes(FRA_young_noise_M.Bayes,FRA_old_noise_M.Bayes,...
      {'Males-Young','Aging'},{'Tones','20','10','0'},...
'Bayes') ;
)

BayesStats.YoungMF =  CompareBayes(FRA_young_noise_M.Bayes,FRA_young_noise_F.Bayes,...
        {'Males','Female'},{'Tones','20','10','0'},...
'Bayes') ;


BayesStats.OldYoungF =  CompareBayes(FRA_young_noise_F.Bayes,FRA_old_noise_F.Bayes,...
        {'Females-Young','Aging'},{'Tones','20','10','0'},...
'Bayes') ;


