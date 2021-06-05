

Old_N = Flatten_Corr(Flatten_Corr(Tonebox_FRA_stats.NCorr));
Old_S = Flatten_Corr(Flatten_Corr(Tonebox_FRA_stats.LCorr));
Young_S =Flatten_Corr(Flatten_Corr(Control_FRA_stats.LCorr));
Young_N = Flatten_Corr(Flatten_Corr(Control_FRA_stats.NCorr));


% output structure
Corrs = []
Corrs.QuietNoise =[]
Corrs.Level = []

% Quiet Noise Comparision
% comparing tone in quiet to same tone in noise
[Corrs.QuietNoise.NCorr.Stats,Corrs.QuietNoise.NCorr.Main] = ...
             CompareOldYoung(Old_N(:,1:2),Young_N(:,1:2),{'99';'20'})
[Corrs.QuietNoise.SCorr.Stats,Corrs.QuietNoise.SCorr.Main] =...
    CompareOldYoung(Old_S(:,1:2),Young_S(:,1:2),{'99';'20'})


% Lowering SNR comparison
% comparing how lowering Tone level affects correlations
[Corrs.Level.NCorr.Stats,Corrs.Level.NCorr.Main] = ...
             CompareOldYoung(Old_N(:,2:4),Young_N(:,2:4),{'20';'10';'0'})
[Corrs.Level.SCorr.Stats,Corrs.Level.SCorr.Main] =...
    CompareOldYoung(Old_S(:,2:4),Young_S(:,2:4),{'20';'10';'0'})

clear Old_N Old_S Young_N Young_S