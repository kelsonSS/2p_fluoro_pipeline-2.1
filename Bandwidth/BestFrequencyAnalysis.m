function [BFs,BF_Field]= BestFrequencyAnalysis(passive)

DF_Levels = passive.df_by_level;
Sig_Levels = passive.df_by_level_sig;

for expt = 1:length(DF_Levels)

freqs = passive.handles{expt}.uFreqs;

Sig_freqs = squeeze(DF_Levels{expt}(1,:,:) .* Sig_Levels{expt}(1,:,:));

[~,BF_expt] =  max(Sig_freqs);
 
BFs{expt} = freqs(BF_expt);

BF_Field(expt) = squeeze(mean(BFs{expt})); 

end 



