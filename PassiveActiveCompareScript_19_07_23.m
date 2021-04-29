




% WARNING: this is a hack that only works with current data set 
% create more permanent solution by looping over data dir for passive file 
% extracting corresponding active file and ensuring the distance between points 
% (pdist2) is close to average dist. could try find peaks of the derivative
% of distnace between points any time there is a change that would indicate
% a mismatch 

size_diff = length(mean_passive) - length(mean_active);

mean_passive = mean_passive(:,1:end-size_diff);

mean_passive = mean_passive';
mean_active = mean_active';


%% sorting and normalizing response 
[passive_max,passive_timing] = max(mean_passive,[],2);
[active_max,active_timing] = max(mean_active,[],2);

mean_active = mean_active./active_max;
mean_passive = mean_passive./passive_max;

[~,passive_order] = sort(passive_timing);
[~,active_order] = sort(active_timing);   


%% base plots

figure;imagesc(mean_passive);title('Average Timecourse passive')
FramesToSeconds

figure;imagesc(mean_active);title('Average Timecourse active')
FramesToSeconds


%% passive ordered plots
mean_passive_p = mean_passive(passive_order,:);
mean_active_p = mean_active(passive_order,:);


figure;imagesc(mean_passive_p);title('Average Timecourse Passive-Passive ordered')
FramesToSeconds

figure;imagesc(mean_active_p);title('Average Timecourse Active-Passive ordered')
FramesToSeconds

%% active ordered plots

mean_passive_a = mean_passive(active_order,:);
mean_active_a = mean_active(active_order,:);


figure;imagesc(mean_passive_a);title('Average Timecourse Passive-Active ordered')
FramesToSeconds


figure;imagesc(mean_active_a);title('Average Timecourse Active-Active ordered')
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











