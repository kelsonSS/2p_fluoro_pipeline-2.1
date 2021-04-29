
Old_Neurons   = squeeze(FRA_old_noise.Bayes.ClassesTonesLossTotal(:,end,:));
Young_Neurons = squeeze(FRA_young_noise.Bayes.ClassesTonesLossTotal(:,end,:));

Old_Timing = squeeze(FRA_old_noise.Bayes.TimeLossTotal(3:end,end,:));
Young_Timing = squeeze(FRA_young_noise.Bayes.TimeLossTotal(3:end,end,:));

Old_SNR = squeeze(FRA_old_noise.Bayes.ClassesTonesLossLvl(80,:,1,:));
Young_SNR =  squeeze(FRA_young_noise.Bayes.ClassesTonesLossLvl(80,:,1,:));

NeuronNumberIDs = createNumberedIDs(100,10);
TimingIDs = createNumberedIDs(150,10);

NN = struct()
[NN.stats, NN.main] =  CompareOldYoung(Old_Neurons',Young_Neurons',NeuronNumberIDs,true);

TT = struct()
[TT.stats, TT.main] = CompareOldYoung(Old_Timing',Young_Timing',TimingIDs,true);

SNR = struct()
[SNR.stats,SNR.main] = CompareOldYoung(Old_SNR,Young_SNR,{'99','20','10','0'})
