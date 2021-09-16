function sig = CorrelationsAnalysisByAnimal(s1,s2,type,GroupNames,Levels)
% CorryAnalysis(S1,S2,type,GroupNames,level)
% Combines Compare2Anova and PlotScatterLevels for Correlations objects
% from CorrelationsByAnimal

title_str = sprintf('%s_%s-%s',GroupNames{1},GroupNames{2},type);  
 

[sig.stats,sig.main_effects] =  Compare2Anova(s1,s2,GroupNames,Levels,type) ; 

close(gcf);
% Plotting only the effects of age for passive paper 

if strmatch(type,'Cell') % if by cell- return here 
 return 
end 

% by animal analysis begins here 
if iscell(s1)
    PlotScatterLevels(GroupNames,"All",title_str,cellfun(@nanmean_all,s1),cellfun(@nanmean_all,s2) )
else  
    PlotScatterLevels(GroupNames,'All',title_str,mean(s1),mean(s2))
end
 
 end
 

function out = nanmean_all(x)

out = nanmean(x(:));
end 





