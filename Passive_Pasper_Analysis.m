function Passive_Paper_Analysis(TN,Savename)

% takes an experiment object and runs all of the analyses used for the
% passive paper 

TN.BD = BandwidthAnalysis(TN,'RFS',1); n


%Passive GC
GCanalTrialsModBalanced_TN(TN,'SNR')
GCanalTrialsModBalanced_TN(TN,'Tones')
GCanalTrialsModBalanced_TN(TN,'Noise')
GCanalTrialsModBalanced_TN(TN,'Offset')
GCanalTrialsModBalanced_TN(TN,'Off')

BayesClassifierPassive(TN,SaveName)


%
