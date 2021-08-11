function x = PlotGroupedErrorBars(Grouped_means, Grouped_CI, calc_flag) 

if exist('calc_flag','var') && calc_flag
    [Grouped_means,Grouped_CI] = Analyze(Grouped_means,Grouped_CI);
end     


h = bar(Grouped_means, 'grouped');
hold on
% Find the number of groups and the number of bars in each group
[ngroups, nbars] = size(Grouped_means);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, Grouped_means(:,i), Grouped_CI(:,i), 'k', 'linestyle', 'none');
end
hold off




function [Grouped_means,Grouped_CI] = Analyze(x1,x2)


Grouped_means = [nanmean(x1);nanmean(x2)];

Grouped_CI = [ nanstd(x1) / sqrt(size(x1,1)) * 1.96 ;  nanstd(x2) / sqrt(size(x2,1)) * 1.96];
      
          
          
          





