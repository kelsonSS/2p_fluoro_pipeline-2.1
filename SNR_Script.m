[SNR.Anova,~,SNR.GCstats]= anova1(SNR.GCnumbers,[],'off');
[SNR.MC_stats] = multcompare(SNR.GCstats,[],'off');

[SNR_Tn.Anova,~,SNR_Tn.GCstats]= anova1(SNR_Tn.GCnumbers,[],'off');
[SNR_Tn.MC_stats] = multcompare(SNR_Tn.GCstats,[],'off');


[SNR_Ns.Anova,~,SNR_Ns.GCstats]= anova1(SNR_Ns.GCnumbers,[],'off');
[SNR_Ns.MC_stats] = multcompare(SNR_Ns.GCstats,[],'off');