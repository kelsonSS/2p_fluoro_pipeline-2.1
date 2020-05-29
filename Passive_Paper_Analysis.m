function Passive_Paper_Analysis(TN,Savename)

% takes an experiment object and runs all of the analyses used for the
% passive paper 
savepath = '\\vault3\data\kelson';

%TN.BD = BandwidthAnalysis(TN,'RFS',1); 

% need to add Correlation Analysis noise correlations, Signal Correlations
TN.Corr = Correlations(TN);

%Passive GC
TN.GC.SNR =GCanalTrialsModBalanced_TN(TN,'SNR');
TN.GC.Tone =GCanalTrialsModBalanced_TN(TN,'Tones');
TN.GC.Noise=GCanalTrialsModBalanced_TN(TN,'Noise');

save(fullfile(savepath,Savename),'TN','-v7.3')
try % these two won't be in quiet expts. 
TN.GC.Ramping = GCanalTrialsModBalanced_TN(TN,'Offset');
TN.GC.NoiseOFF = NGCanalTrialsModBalanced_TN(TN,'Off');
catch
end 


save(fullfile(savepath,Savename),'TN','-v7.3')
TN.Corrs = Corrs_by_temporal_timing(TN);
save(fullfile(savepath,Savename),'TN','-v7.3')

% Bayes Figures and File Saving 
BayesClassifierPassive(TN,Savename)


%
