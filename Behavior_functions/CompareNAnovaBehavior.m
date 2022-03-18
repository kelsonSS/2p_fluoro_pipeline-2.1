function [all_stats] = Compare2AnovaBehavior(behave,GroupNames,levels,SaveName,varnames)

% extends Compare2Anova to work with arbitrarily many groups and Levels

% 
% behave is now a group x level cell with each cell containing a vector
% corresponding to the target variable at the correct group and level

if ~ exist('SaveName','var')
    SaveName = [];
end 

if ~exist('varnames','var')
    varnames = {'Group','Performance'};
end 


b_flat = cellfun(@(x) x(:), behave,'UniformOutput',false);
  
group_idx = getIndex(behave,1,GroupNames);
lvl_idx = getIndex(behave,2,levels);


all_data = vertcat(b_flat{:});




% create final indices

[x,main_effects,stats,z]= anovan(all_data,{group_idx,lvl_idx},'model','interaction','varnames',varnames);



[means,CIs] = getMeansAndCIs(behave);

figure
PlotGroupedErrorBars(means,CIs)
xticklabels(GroupNames)
title(SaveName,'Interpreter','none')

% if SaveName 
%     suptitle(sprintf('Main Effects for: %s', SaveName'))
%     saveas(gcf,sprintf('%s -MainEffects.pdf',SaveName) )
% end 


main_effects_P_value = main_effects(2:3,end);
 my_maineffectsplot(all_data,{group_idx,lvl_idx},...
                    main_effects_P_value,...
                    {'age','level'})

                if SaveName 
     suptitle(sprintf('Main Effects for: %s', SaveName'))
     saveas(gcf,sprintf('%s -MainEffectsPlot.pdf',SaveName) )
 end 


if length(levels) > 1
    figure
    all_stats.comparisions = multcompare(stats,'Dimension',[1 2]);
end

all_stats.stats = stats;
all_stats.main_effect = main_effects;


function full_idx = getIndex(data,dim,IDs)

 
    if dim == 1
        full_idx = BuildIndex(data',IDs);
        full_idx = full_idx'
       
    else
        full_idx = BuildIndex(data,IDs);
    end 
      full_idx = vertcat(full_idx{:});
    
  
    
    function full_idx = BuildIndex(data,IDs)       
        
    n_groups = size(data,1);
    n_levels = size(data,2);
    
    full_idx = cell(n_groups,n_levels);
    for row = 1:n_groups
        for col = 1:n_levels
            n_expts = length(data{row,col});
    
            full_idx(row,col) = {repmat(IDs(col),n_expts,1)};
        end 
    end  
 
        
    





function [Means,CIs] =  getMeansAndCIs( b1)
        
        Means = cellfun(@nanmean,b1);
        CIs =  cellfun(@getCI_95,b1);
            


function CI = getCI_95(data)
      
        CI = nanstd(data,[],2)./sqrt(size(data,2) * 1.96);              


