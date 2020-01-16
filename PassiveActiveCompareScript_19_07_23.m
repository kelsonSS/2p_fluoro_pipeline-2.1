

DFF_passive = BehaviorPassive.DFF;

passive_freq_idx = BehaviorPassive.FreqLevelOrder{:,1} == 11314;

passive_expt_idx = BehaviorPassive.Experiment_list == 1 |...
                   BehaviorPassive.Experiment_list == 2;

DFF_passive = DFF_passive(:,:,passive_expt_idx);

DFF_passive  = squeeze(nanmean(DFF_passive,2));



DFF_active = []
for ii =1:2:3
DFF_temp = squeeze(nanmean(test{ii}.DFF,2));
DFF_active = cat(2,DFF_active,DFF_temp);
end
DFF_active = DFF_active(31:120,:);


% WARNING: this is a hack that only works with current data set 
% create more permanent solution by looping over data dir for passive file 
% extracting corresponding active file and ensuring the distance between points 
% (pdist2) is close to average dist. could try find peaks of the derivative
% of distnace between points any time there is a change that would indicate
% a mismatch 

size_diff = length(DFF_passive) - length(DFF_active);

DFF_passive = DFF_passive(:,1:end-size_diff);

DFF_passive = DFF_passive';
DFF_active = DFF_active';


%% sorting and normalizing response 
[passive_max,passive_timing] = max(DFF_passive,[],2);
[active_max,active_timing] = max(DFF_active,[],2);

DFF_active = DFF_active./active_max;
DFF_passive = DFF_passive./passive_max;

[~,passive_order] = sort(passive_timing);
[~,active_order] = sort(active_timing);   


%% base plots

figure;imagesc(DFF_passive);title('Average Timecourse passive')
FramesToSeconds

figure;imagesc(DFF_active);title('Average Timecourse active')
FramesToSeconds


%% passive ordered plots
DFF_passive_p = DFF_passive(passive_order,:);
DFF_active_p = DFF_active(passive_order,:);


figure;imagesc(DFF_passive_p);title('Average Timecourse Passive-Passive ordered')
FramesToSeconds

figure;imagesc(DFF_active_p);title('Average Timecourse Active-Passive ordered')
FramesToSeconds

%% active ordered plots

DFF_passive_a = DFF_passive(active_order,:);
DFF_active_a = DFF_active(active_order,:);


figure;imagesc(DFF_passive_a);title('Average Timecourse Passive-Active ordered')
FramesToSeconds


figure;imagesc(DFF_active_a);title('Average Timecourse Active-Active ordered')
FramesToSeconds

figure;histogram(active_timing - passive_timing);
xticks([-60,-30,0, 30,60])
xticklabels({'-2','-1','0','1','2'})
xlabel('Difference in Seconds (active - passive) ')   
ylabel('count')
title('Difference between Active and Passive Peak')

%% difference plots 

figure;scatter(passive_timing(passive_order),1:length(passive_order))
hold on 
scatter(active_timing(passive_order),1:length(passive_order))
axis ij


figure;scatter(active_timing(active_order),1:length(active_order))
hold on 
scatter(passive_timing(active_order),1:length(active_order))
axis ij











