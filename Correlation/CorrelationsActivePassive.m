function out = CorrelationsActivePassive(TNBehavior,SaveName,BehaviorType,PlotColor)
if ~exist('PlotColor','var')
    PlotColor = 'k'
end 

if ~exist('SaveName','var')
    SaveName = [];
end 

if ~exist('BehaviorType','var')
    BehaviorType = 'Hit';
end 


[active,passive,TNBehavior] = getActivePassiveCorrs(TNBehavior,BehaviorType);




close_idx_all  = CompareBestFrequencyToBehaviorTarget(TNBehavior);

close_idx = cell2mat(close_idx_all);


out =  FullCompareCorrs(passive,active,...
                             close_idx,PlotColor,...
                           [SaveName '_' BehaviorType]);
                       
                       
                       