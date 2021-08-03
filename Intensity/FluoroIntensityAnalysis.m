function sig = FluoroIntensityAnalysis(s1,s2,type,GroupNames,Levels)
% CorryAnalysis(S1,S2,type,GroupNames,level)
% Combines Compare2Anova and PlotScatterLevels for Correlations objects
% from CorrelationsByAnimal

  s1_cell = PlotFluoroCDF(s1,type,Levels,GroupNames{1});
  s2_cell = PlotFluoroCDF(s2,type,Levels,GroupNames{2});

  [~,sig.cell] =ttest2(s1_cell,s2_cell);
 
 
 %byAnimal
     
  s1_animal = PlotFluoroByAnimal(s1,type);

  s2_animal = PlotFluoroByAnimal(s2,type);

  
 [~,sig.animal] = ttest2(s1_animal,s2_animal);

title_str = sprintf('%s_%s-DFTone-%s%s',GroupNames{1},GroupNames{2},type,num2str(Levels))  
 

 PlotScatterLevels(GroupNames,num2str(Levels),title_str,s1_animal,s2_animal)
 





