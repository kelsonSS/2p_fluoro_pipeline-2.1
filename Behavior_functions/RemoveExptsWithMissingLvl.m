function [Behavior,BehaviorSheet] = RemoveExptsWithMissingLvl(Behavior,BehaviorSheet,lvls_to_use)

% this function is part of TNDetectionAnalysis and is used to only select
% from experiments that have only at most 1 level missing from experiments

n_levels = length(lvls_to_use);
good_idx =cellfun(@(x) length(intersect(x,lvls_to_use)),BehaviorSheet.SNRs);
good_idx = good_idx >= n_levels-1;


Behavior = Behavior(good_idx,:);
BehaviorSheet = BehaviorSheet(good_idx,:);
