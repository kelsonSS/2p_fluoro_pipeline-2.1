function sig = CorrelationsAnalysisByCell(s1,s2,GroupNames,Levels,type)

[sig.stats,sig.main_effects] =  Compare2Anova(s1,s2,GroupNames,Levels,type) ; 
