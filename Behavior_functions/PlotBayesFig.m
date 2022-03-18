
function stats = PlotBayesFig(Behavior,field,SaveName)

stats = [];
if ~exist('SaveName','var')
    SaveName = [];
end 

colors = {'b','r','k','g','r'}

BehaveNames = fieldnames(Behavior);
n_groups = length(BehaveNames);
figure
set(gcf,'Position',[ 200,100,500, 800])
hold on 
data = {};
for group_name_idx = 1:n_groups
% extract
    data{group_name_idx} = Behavior.(BehaveNames{group_name_idx}).(field);
     subplot(n_groups,1,group_name_idx)
     plotData(data{group_name_idx},field,colors{group_name_idx})

    title( BehaveNames{group_name_idx})
   
% beautify
    hold on 
    ax = gca;
    plot([ax.XLim(1) ax.XLim(2)], [0 0],'k--') 
    %plot( [30 30],[ax.YLim(1) ax.YLim(2)],'b--') 
    plot([60 60],[ax.YLim(1) ax.YLim(2)],'b--') 
   
    ylabel('Normalized % Correct')
    ylim([-10 23])

%     % comment out for BayesNumbers 
%      xticks(0:30:90)
%     xticklabels(0:1:3)
%   if group_name_idx == n_groups
%     xlabel("Time (S)")
%     xticks(0:30:90)
%     xticklabels(0:1:3)
%   end 

end 



suptitle(strrep(field,'_','-'))

if SaveName
    saveas(gcf,[SaveName '.pdf']);
end 

if contains(field,'Level')
    figure
    [Means, CIs,raw] = MungeData(data);
   
    Levels = {'0 dB','10dB','20dB'};
    PlotGroupedErrorBars(Means',CIs')
    legend(BehaveNames)
    xticklabels(Levels)
    ylabel('% Correct')
    saveas(gcf,[SaveName '-Bars.pdf'])
    
  stats =  CompareNAnovaBehavior(raw,BehaveNames,Levels,SaveName);
    
end


function [means, CIs, expt_mu ] = MungeData(data)


means = [];
CIs = [];
exp_mu = {};
for group_idx = 1:length(data)
    
    [means(group_idx,:),CIs(group_idx,:),expt_mu(group_idx,:) ] = getRelData(data{group_idx});
    
    
end


function [means,CIs,raw] = getRelData(data)
% data is a level x time x experiment x rep matrix
% output is level x experiment matrix 

% find max response, find baseline 
n_levels = size(data,1);

means = zeros(n_levels,1);
CIs = zeros(n_levels,1);

window_start = 5;
window_end = 60;
raw = {};

for lvl = 1 :n_levels 
% baseline correct
    corrected_data = baseline_corrected_average(data(lvl,2:end,:,:));
  
 % find best performance during tone period !!!warning!!! currently hardcoded    
    
    [~,best_idx ] = max( nanmean(corrected_data(window_start:window_end,:),2) );
      % adjust index and get get mean and CI from corrected index 
    
   
    best_idx = best_idx + window_start -1;
    
    raw{lvl} = corrected_data(best_idx,:);
    means(lvl) = nanmean(corrected_data(best_idx,:));
    CIs(lvl) = getCI(corrected_data(best_idx,:));
    
    
end


function CI = getCI(rel_max_data)
n = length(rel_max_data);

CI = squeeze(nanstd(rel_max_data))/ sqrt(n) * 1.96 ;











function plotData(data,field,color)
 
colors = {'r','g','b'}

if ndims(data) == 3
    data_shape = size(data);
    data = reshape(data,1,data_shape(1),data_shape(2),data_shape(3));
    n_levels = 1 ;
    colors = {color};
    
elseif ndims(data) == 4
    n_levels = size(data,1);
else 
    error('%s is not a valid field to plot',field)
end
    
    
hold on
for lvl = 1 :n_levels 
% baseline correct
    corrected_data = baseline_corrected_average(data(lvl,2:end,:,:));
% plot     
   
    plotShadedErrorBar(corrected_data,colors{lvl} )
    
end



function corrected_mu = baseline_corrected_average(d)

baseline_frames = 30;
trialdur = size(d,1);

mu = squeeze(reshape(d,1,size(d,2),[]));

mu_baseline = repmat( nanmean(mu(1:baseline_frames,:)),[trialdur,1]);
corrected_mu = (mu - mu_baseline) * 100;



