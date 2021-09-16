function CompareCorrsByDistance(C1,C2,SaveName) 



PlotCorrs(C1.LCorr,C2.LCorr,C1.BinSize,[SaveName '-', 'SignalCorr'])
PlotCorrs(C1.NCorr,C2.NCorr,C1.BinSize,[SaveName '-', 'NoiseCorr'])




function PlotCorrs(C1,C2,BinSizeMicrons,SaveName)



x = [0:BinSizeMicrons:(size(C1.Cell_mu,1)-1)*BinSizeMicrons];
levelTitles = {'Tones','+20','+10','0'}

figure
for ii =1:4
subplot(2,2,ii)
shadedErrorBar(x,C1.Cell_mu(:,ii),C1.Cell_CI(:,ii),'b'); hold on
shadedErrorBar(x,C2.Cell_mu(:,ii),C2.Cell_CI(:,ii),'k'); 
title(levelTitles{ii})
ylim([0 .5])
ylabel('Correlation (AU)')
xlabel('Distance (um)')
set(gcf,'Position',[660 450 560 520])
    
end 
suptitle(SaveName)

saveas(gcf,sprintf('%s-ByDistance.pdf', SaveName))


    


