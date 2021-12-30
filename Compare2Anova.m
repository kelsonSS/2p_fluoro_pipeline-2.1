function [all_stats,main_effects] = Compare2Anova(g1,g2,GroupNames,levels,SaveName,varnames,full_levels)

%  Compare2Anova(g1,g2,GroupNames,levels,SaveName,varnames,full_levels)

% this function will create an anova to compare two groups
% at different levels 
%
% input 
% old/young  matrices in a samples X levels format
% levels - 1 x N cell of different levels used {'90';'20';'10;'0'}

if ~exist('varnames','var')
    varnames = {'age','level'};
end 

if ~exist('full_levels','var')
    full_levels = false;
end 

if ~exist('SaveName','var')
    SaveName = [];
else 
    SaveName = [GroupNames{1}, GroupNames{2},'-',SaveName];
end 


% plot error-bars 
if iscell(g1) || iscell (g2)
   g1 = cell2mat(cellfun(@(x) x(:), g1,'UniformOutput',false)')';
   g2 = cell2mat(cellfun(@(x) x(:), g2,'UniformOutput',false)')';
    
end
    
Means =[nanmean(g1,2),nanmean(g2,2)]
CIs = [ nanstd(g1,[],2)./sqrt(size(g1,2)) * 1.96, ...
        nanstd(g2,[],2)./sqrt(size(g2,2)) * 1.96  ];   

figure;PlotGroupedErrorBars(Means , CIs )
legend(GroupNames)
xticklabels(levels)
title(SaveName,'interpreter','none')
xlabel('Levels')
ylabel('A.U')
saveas(gcf,sprintf('%s_ComparisonBars.pdf',SaveName))

% create indexing for Anova

g1_idx = repmat({GroupNames{1}}, numel(g1),1);

g2_idx = repmat({GroupNames{2}},numel(g2),1);
age_idx = cat(1,g1_idx,g2_idx);

g1_lvl_idx = repmat(levels',1,size(g1,2));

g1 = g1(:);
g2 = g2(:);

all = cat(1,g1,g2);

if ~ full_levels
lvl_idx = repmat(levels,1,length(all)/length(levels) ) ;
lvl_idx = lvl_idx(:);
else 
    lvl_idx = full_levels;
end 


clean_idx = ~isnan(all);
all = all(clean_idx);
lvl_idx = lvl_idx(clean_idx);
age_idx = age_idx(clean_idx);

[x,main_effects,stats,z]= anovan(all,{age_idx,lvl_idx},'model','interaction','varnames',varnames);



% if SaveName 
%     suptitle(sprintf('Main Effects for: %s', SaveName'))
%     saveas(gcf,sprintf('%s -MainEffects.pdf',SaveName) )
% end 


main_effects_P_value = main_effects(2:3,end);
 my_maineffectsplot(all,{age_idx,lvl_idx},...
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
end 
% % 
% % 
% % if SaveName
% %     saveas(gcf,sprintf('%s -Comparision.pdf',SaveName) )
% % end 


