function sig = CorrelationsAnalysisByAnimal(s1,s2,type,GroupNames,Levels)
% CorryAnalysis(S1,S2,type,GroupNames,level)
% Combines Compare2Anova and PlotScatterLevels for Correlations objects
% from CorrelationsByAnimal

title_str = sprintf('%s_%s-Corr-%s',GroupNames{1},GroupNames{2},type);  
 
[sig.stats,sig.main_effects] =  Compare2Anova(s1,s2,GroupNames,Levels,title_str) ; 

 PlotScatterLevels(GroupNames,Levels,title_str,s1,s2)
 
 
 end 





