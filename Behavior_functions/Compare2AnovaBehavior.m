function [all_stats,main_effects] = Compare2AnovaBehavior(behave1,behave2,GroupNames,levels,SaveName,varnames)

% extends Compare2Anova to work with behavioraldata

% behave1 and 2 are now a N x 1 cell with each cell containing a vector
% corresponding to the target variable at the Nth 'level'

if ~ exist('SaveName','var')
    SaveName = [];
end 

if ~exist('varnames','var')
    varnames = {'Age','Performance'};
end 


b1_flat = cell2mat(cellfun(@(x) x(:), behave1,'UniformOutput',false)');
  
b1_levels = getLevels(behave1,levels);
b1_group =  repmat(GroupNames(1),length(b1_flat),1);

b2_flat = cell2mat(cellfun(@(x) x(:), behave2,'UniformOutput',false)');

b2_levels = getLevels(behave2,levels);

b2_group =  repmat(GroupNames(2),length(b2_flat),1);

% create final indices
all = cat(1,b1_flat,b2_flat);
group_idx = cat(1,b1_group,b2_group);
lvl_idx = cat(1,b1_levels,b2_levels);

[x,main_effects,stats,z]= anovan(all,{group_idx,lvl_idx},'model','interaction','varnames',varnames);



% if SaveName 
%     suptitle(sprintf('Main Effects for: %s', SaveName'))
%     saveas(gcf,sprintf('%s -MainEffects.pdf',SaveName) )
% end 


main_effects_P_value = main_effects(2:3,end);
 my_maineffectsplot(all,{group_idx,lvl_idx},...
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


function full_lvls = getLevels(b,levels)


lvl_idx = [];

for curr_lvl = 1:length(b)
    
  curr_idx =  repmat(curr_lvl,length(b{curr_lvl}),1);
    
  
 lvl_idx = cat(1,lvl_idx,curr_idx);
end 

full_lvls = levels(lvl_idx)';
    
    


