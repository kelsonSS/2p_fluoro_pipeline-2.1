% %JKM Processing Script
% %Get Data Dirs
% dataDir=ans.DataDirs;
% 
% %%
% %%Fluoro_to_Table
% Fluoro_to_Table(dataDir)
% DF=ans;
%%
%%Get Clusters
plotAllDFFs(DF)
Cluster_DF(DF)

%%
%%Get Bandwidth
BandwidthAnalysis(DF, 'Significant',0,'Pos',.5)
BandwidthData=ans
% oldBandwidthAnalysis(DF, 'Significant',0,'Pos',.5)
% oldBandwidthAnalysis(DF, 'Significant',0,'Neg',.5)

%%
%%Get Cell Response Timing
% PlotCellResponseTiming(DF)
TemporalAnalysis(DF)
FluoroResponseTimes=ans

%% Plot Mean & Max Flouresence
PlotFluoroCDF(DF,'max')
MaxFluoro=ans
PlotFluoroCDF(DF,'mean')
MeanFluoro=ans

%% Plot Signal & Noise Correlations
Correlations(DF)
CorrelationStats=ans