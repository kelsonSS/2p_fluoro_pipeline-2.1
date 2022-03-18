function histogramLevels(data,edges,levelNames,SaveName)

if ~exist('SaveName','var')
    SaveName = []'
end 

n_levels = length(levelNames);

figure

for lvl = 1:n_levels
    
    subplot(n_levels,1,lvl)
    c = histcounts(abs(data{lvl}),edges) / length(data{lvl}) * 100
    plot(edges(2:end),cumsum(c), 'LineWidth',2);
    title(levelNames{lvl})
    set(gca,'xscale','log')
    ylim([0 100])
    ylabel('%')
end 


if SaveName
    saveas(gcf,[SaveName '-CDF_Levels.pdf'])
end 


