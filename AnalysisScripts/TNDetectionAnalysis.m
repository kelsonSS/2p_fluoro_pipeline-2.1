% BehaviorAnalysis Script

BehaviorList = FindPairedBehaviorExperiments;
TNBehavior = Fluoro_to_Table_Behavior(BehaviorList);

% remove any rows that didn't properly get extracted 
has_behavior_idx = ~cellfun(@isempty,TNBehavior(:,2));
TNBehavior = TNBehavior(has_behavior_idx,:);

% Behavior has some files with different amount of PreStim Silence
% trim them all to all be consistent with 1 second prestim silence 
TNBehavior(:,2) =  TNBehavior_Fix_Timings(TNBehavior(:,2));

old_idx = contains(lower(TNBehavior(:,3)),'ia');

TNBehaviorYoung = TNBehavior(~old_idx,:);
TNBehaviorOld = TNBehavior(old_idx,:);

%% create output struct
BehaviorResults = struct();

BehaviorResults.Young.CellInfo = CreateBehaviorCellList(TNBehaviorYoung);
BehaviorResults.Old.CellInfo = CreateBehaviorCellList(TNBehaviorOld);

AnimalIDS = cellfun(@(x) strsplit(x,filesep),TNBehavior(:,3),'UniformOutput',0);
 AnimalIDS = unique( cellfun(@(x) x{4},AnimalIDS,'UniformOutput',0));
 BehaviorResults.AnimalID = AnimalIDS;
 clear AnimalIDS
 
%% clustering 
BehaviorResults.Young.Clust.Passive = Cluster_DF(TNBehaviorYoung(:,1),'K-means');
BehaviorResults.Young.Clust.Active = Cluster_DF(TNBehaviorYoung(:,2),'K-means',9,...
'normalized',BehaviorResults.Young.Clust.Passive.Centroids);

BehaviorResults.Old.Clust.Passive = Cluster_DF(TNBehaviorOld(:,1),'K-means');
BehaviorResults.Old.Clust.Active = Cluster_DF(TNBehaviorOld(:,2),'K-means',9,...
'normalized',BehaviorResults.Old.Clust.Passive.Centroids);

BehaviorResults.Combined.Clust.Passive = Cluster_DF(TNBehavior(:,1),'K-means',10);
BehaviorResults.Combined.Clust.Active = Cluster_DF(TNBehavior(:,2),'K-means',6,...
'normalized',BehaviorResults.Combined.Clust.Passive.Centroids);


%% plot Cluster transitions 
% Young_passive_clusters = BehaviorResults.Young.Clust.Passive.Clusters(...
%      BehaviorResults.Young.Clust.Passive.Neuron);
%  Young_active_clusters = BehaviorResults.Young.Clust.Active.Clusters(...
%      BehaviorResults.Young.Clust.Active.Neuron);
%  
%  Old_passive_clusters = BehaviorResults.Old.Clust.Passive.Clusters(...
%      BehaviorResults.Old.Clust.Passive.Neuron);
%  Old_active_clusters = BehaviorResults.Old.Clust.Active.Clusters(...
%      BehaviorResults.Old.Clust.Active.Neuron);
% 
%  % plot Young
%  trans_mat_young = histcounts2(Young_passive_clusters,Young_active_clusters(1:2892));
% trans_mat_young = trans_mat_young ./ sum(trans_mat_young(:)) ;
%  figure; imagesc(trans_mat_young)
 % Plot Old 
 
 %%
 
  behavior_path ='Z:\Kelson\TNDetectionAnalysis';
% load in behavior 
for ii = 1:length(BehaviorResults.AnimalID)
    curr_id = BehaviorResults.AnimalID{ii};
    
    if contains(curr_id,'IA','IgnoreCase',true)
    BehaviorResults.Old.Behavior.( sprintf('animal_%s',curr_id) ) = ...
        load(fullfile(behavior_path,curr_id));
    else 
         BehaviorResults.Young.Behavior.( sprintf('animal_%s',curr_id) ) = ...
        load(fullfile(behavior_path,curr_id));
    end 
end 

TN_passive =Munge_DF(TNBehavior(:,1));
TN_active = Munge_DF(TNBehavior(:,2));

% analyze behavior
BehaviorResults.Old.Behavior.group = MungeBehaviorGroupData(BehaviorResults.Old.Behavior);
BehaviorResults.Young.Behavior.group = MungeBehaviorGroupData(BehaviorResults.Young.Behavior);

%Plot group Behavior 
PlotGroupedData(BehaviorResults.Old.Behavior.group)
PlotGroupedData(BehaviorResults.Young.Behavior.group)


% analyse Passive Data

