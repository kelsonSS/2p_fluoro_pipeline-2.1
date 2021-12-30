
function PlotBayesFig(Behavior,field,SaveName)

if ~exist('SaveName','var')
    SaveName = [];
end 

colors = {'k','g','r','b','c','m','y'}

BehaveNames = fieldnames(Behavior);
n_groups = length(BehaveNames)
figure
set(gcf,'Position',[ 200,100,500, 800])
hold on 
for group_name_idx = 1:n_groups
% extract
    data = Behavior.(BehaveNames{group_name_idx}).(field);
% baseline correct
    corrected_data = baseline_corrected_average(data);
% plot     
    subplot(n_groups,1,group_name_idx)
    plotShadedErrorBar(corrected_data,colors{group_name_idx} )
    title( BehaveNames{group_name_idx})
   
% beautify
    hold on 
    ax = gca;
    plot([ax.XLim(1) ax.XLim(2)], [0 0],'k--') 
    %plot( [30 30],[ax.YLim(1) ax.YLim(2)],'b--') 
    plot([60 60],[ax.YLim(1) ax.YLim(2)],'b--') 
   
    ylabel('Normalized % Correct')
    ylim([-10 23])

    % comment out for BayesNumbers 
     xticks(0:30:90)
    xticklabels(0:1:3)
  if group_name_idx == n_groups
    xlabel("Time (S)")
    xticks(0:30:90)
    xticklabels(0:1:3)
  end 

end 



suptitle(strrep(field,'_','-'))

if SaveName
    saveas(gcf,[SaveName '.pdf']);
end 

function corrected_mu = baseline_corrected_average(d)

baseline_frames = 4;
trialdur = size(d,1);

mu = squeeze(nanmean(d,3));

mu_baseline = repmat( nanmean(mu(1:baseline_frames,:)),[trialdur,1]);
corrected_mu = (mu - mu_baseline) * 100;



