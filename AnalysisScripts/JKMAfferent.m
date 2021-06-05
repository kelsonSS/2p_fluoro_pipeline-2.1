%JKM Processing Script
%Get Data Dirs
dataDir=DF.DataDirs;

%%
%%Fluoro_to_Table
Fluoro_to_Table(dataDir)
DF=ans;
%%
%%Get Clusters
plotAllDFFs(DF)
Cluster_DF(DF)

%%
%Get Bandwidth
BandwidthAnalysis(DF, 'Significant',0,'Pos',.5)
BandwidthData=ans

%%
% Get BF Distribution
PlotBFDistribution(DF)
BFDistribution=ans

%%
%Get Cell Response Timing
% PlotCellResponseTiming(DF)
TemporalAnalysisByAnimal(DF)
FluoroResponseTimes=ans

%% Plot Mean & Max Flouresence
PlotFluoroCDF(DF,'max')
MaxFluoro=ans
PlotFluoroCDF(DF,'mean')
MeanFluoro=ans

%% Plot Signal & Noise Correlations
Correlations(DF)
CorrelationStats=ans