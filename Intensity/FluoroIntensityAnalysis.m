function sig = FluoroIntensityAnalysis(s1,s2,type,GroupNames,Levels)
% CorryAnalysis(S1,S2,type,GroupNames,level)
% Combines Compare2Anova and PlotScatterLevels for Correlations objects
% from CorrelationsByAnimal
if ~exist('Levels','var')
    Levels = [];
end 
GroupNames = cellfun(@(x) strrep(x,'-','_'),GroupNames,'UniformOutput',0);


SaveName = sprintf('%s_%s',GroupNames{1},GroupNames{2}); 

  s1_cell = PlotFluoroCDF(s1,type,Levels); close(gcf);
  s2_cell = PlotFluoroCDF(s2,type,Levels); close(gcf);

 figure
cdfplot(s1_cell)
hold on 
cdfplot(s2_cell)
legend(GroupNames{1},GroupNames{2})
xlabel('Fluorescence DF/F')
ylabel('Culmulative Fraction')

title_str = sprintf('FluoroCDF-%s-%s-Cell',type,num2str(Levels))
    title(title_str)
  
        saveas(gcf, sprintf('%s-%s.pdf', SaveName, title_str) )
  
  
  [~,sig.cell,~,sig.tstatCell] =ttest2(s1_cell,s2_cell);
 
 
 %byAnimal
     
  s1_animal = PlotFluoroByAnimal(s1,type); close(gcf)

  s2_animal = PlotFluoroByAnimal(s2,type); close(gcf)

  
 [~,sig.animal,~,sig.tstatAnimal] = ttest2(s1_animal,s2_animal);

title_str = sprintf('%s-FluoroCDF-%s-%s-Animal',SaveName,type,num2str(Levels))  
 

 PlotScatterLevels(GroupNames,num2str(Levels),title_str,s1_animal,s2_animal)
 
% Get Means

 sig.(sprintf('%s_cell_mu', GroupNames{1})) = mean(s1_cell);
 sig.(sprintf('%s_cell_mu', GroupNames{2})) = mean(s2_cell);

 sig.(sprintf('%s_animal_mu', GroupNames{1})) = mean(s1_animal);
 sig.(sprintf('%s_animal_mu', GroupNames{2})) = mean(s2_animal);



