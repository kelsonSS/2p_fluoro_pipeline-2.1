function out= CompareCorrsByDistanceBehavior(Young,Old,SaveName) 

out.Young = Young;
out.Old = Old;

%PlotCorrs(Young.LCorr,Old.LCorr,Young.BinSize,[SaveName '-', 'SignalCorr'])
PlotCorrs(Young.NCorr,Old.NCorr,Young.BinSize,[SaveName '-', 'NoiseCorr'])




function PlotCorrs(C1,C2,BinSizeMicrons,SaveName)



x = [0:BinSizeMicrons:(size(C1.Cell_mu,1)-1)*BinSizeMicrons];
levelTitles = {'+20','+10','0'}

figure
for ii =1:3
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


    


