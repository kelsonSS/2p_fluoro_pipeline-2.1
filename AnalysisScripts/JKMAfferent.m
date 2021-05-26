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
% oldBandwidthAnalysis(DF, 'Significant',0,'Pos',.5)
% oldBandwidthAnalysis(DF, 'Significant',0,'Neg',.5)

%%
%%Get Cell Response Timing
PlotCellResponseTiming(DF)
TemporalAnalysis(DF)

%% Plot Mean & Max Flouresence
PlotFluoroCDF(DF,'max')
PlotFluoroCDF(DF,'mean')

%% Plot Signal & Noise Correlations
Correlations(DF)