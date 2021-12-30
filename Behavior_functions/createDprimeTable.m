function dtable = createDprimeTable(resultsTable)
snrs = [ -10 -5 0 5 10 20]; 

uSNRS = length(snrs);
N_expts = size(resultsTable,1);

dtable = nan(N_expts+1,uSNRS+1);


for expt_idx = 1:N_expts
    for SNR_idx = 1:uSNRS
    
    curr_SNR_id = find(resultsTable.SNRs{expt_idx} == snrs(SNR_idx));

    if curr_SNR_id
   
        dtable(expt_idx,SNR_idx) = resultsTable.dPrimeLevel{expt_idx}(curr_SNR_id);
   end 
    end 
end 

dtable(dtable == Inf) = 3;
dtable(dtable == -Inf) = -3;

plotDprimeFig(dtable,snrs,resultsTable.AnimalID)
saveas(gcf,'AnimalDprimeBySNR.pdf')

plotDprimeFig(dtable(:,[3,5,6,end]),snrs([3,5,6]),resultsTable.AnimalID)
saveas(gcf,'AnimalDprimeBySNR-0_10_20dB.pdf')
PlotCI(dtable(:,[3,5,6]),{'0','10','20'})
xlabel('dB SNR')
ylabel('dPrime')
title('Detection accuracy by level')
saveas(gcf,'AverageDprimeBySNR.pdf')
dtable = dtable(1:end-1,1:end-1);

function plotDprimeFig(dtable,snrs,AnimalID)
uSNRS = length(snrs);
N_expts = size(AnimalID);
f= figure
f.Position = [680   291   817   687];
pcolor(dtable)
xticks(1.5:uSNRS+1.5)
xticklabels(snrs);
yticks(1.5:N_expts+1.5)
yticklabels(AnimalID)
ylabel('Animal')
xlabel('SNR')
title('Dprime By SNR')
colorbar
caxis([0 3])


    
