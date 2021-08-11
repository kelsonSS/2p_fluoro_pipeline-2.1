function my_maineffectsplot(all,groups,p_values,varnames)

% this function plots takes the inputs in the same style that would be put
% into AnovaN and plots the main effects. currently only 2way anova is
% supported

add_scatterplot = false;

h = figure
subplot(1,2,1)
hold on
[Means1,CI1,LevelNames1] = CalculateMeanCI(all,groups{1});

x_vals=PlotGroupedErrorBars(Means1,CI1);
if add_scatterplot
    PlotGroupedScatter(all,groups{1},x_vals)
end 
CleanPlot(varnames{1},LevelNames1,p_values{1})


subplot(1,2,2)
hold on
[Means2,CI2,LevelNames2] = CalculateMeanCI(all,groups{2});

x_vals2 =PlotGroupedErrorBars(Means2,CI2);
if add_scatterplot
    PlotGroupedScatter(all,groups{2},x_vals2)
end 


CleanPlot(varnames{2},LevelNames2,p_values{2})


function [all_mean,all_CI,UIDs] = CalculateMeanCI(allData,IDS)
% find unique groups
UIDs = unique(IDS);
N_IDS = length(UIDs);

% preallocate outputs
all_mean = zeros(N_IDS,1);
all_CI = zeros(N_IDS,1);

for ID_idx = 1:N_IDS
%  calc mean and CI for each group
    group_data = allData(strcmp(UIDs{ID_idx},IDS)); 

    all_mean(ID_idx) = nanmean(group_data);
    all_CI(ID_idx) = CalculateConfidenceInterval(group_data);
    

end 
  

function CleanPlot(title_str,xlabels,p)

title(sprintf('Main Effect of %s:P = %.2d',title_str,p))

xticks(1:numel(xlabels))
xticklabels( xlabels)



function CI = CalculateConfidenceInterval(Data)
    % standard formula to calculate 95% CI 

CI = nanstd(Data) / sqrt( size(Data,1)) *  1.96; 


function PlotGroupedScatter(Data,group_ID,x_vals)

 UID= unique(group_ID);
 IDs = cellfun(@(x) find(strcmp(UID,x)),group_ID);
 
if exist('x_vals','var')
  IDs = x_vals(IDs);
end 
hold on
gscatter(IDs,Data,group_ID)



