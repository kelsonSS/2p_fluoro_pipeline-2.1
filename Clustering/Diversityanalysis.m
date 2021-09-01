function stats =  Diversityanalysis(g1,g2,GroupNames)

g1_animal = ClusterDiversityByAnimal(g1);
g2_animal = ClusterDiversityByAnimal(g2);


scatter_title = sprintf('%s%s-ClusterDiversity-Animal-Scatter',GroupNames{1},GroupNames{2});

PlotScatterLevels(GroupNames,'Group',scatter_title,g1_animal,g2_animal)

[~,stats.p,~,stats.tstat] = ttest2(g1_animal,g2_animal);

figure
PlotGroupedErrorBars(g1_animal,g2_animal,1)
xticklabels(GroupNames)
ylabel('# Active clusters')
title('Cluster-Diversity')

saveas(gcf, sprintf('%s%s-clusterDiversityAnimal.pdf',GroupNames{1},GroupNames{2}) );


g1_cell = Plot_ClusterDiversity(g1.Combined_Classes);close(gcf)
g2_cell = Plot_ClusterDiversity(g2.Combined_Classes);close(gcf)

figure
bar([g1_cell;g2_cell]')
legend(GroupNames)
xlabel('Cluster ID')
ylabel('Fraction of total neurons')
title('Cluster Diversity')
saveas(gcf, sprintf('%s%s-clusterDiversityCell.pdf',GroupNames{1},GroupNames{2}) );





