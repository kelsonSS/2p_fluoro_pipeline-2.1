function [all_stats,main_effects] = Compare2Anova(g1,g2,GroupNames,levels,SaveName,varnames,full_levels)

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
end 


g1_idx = repmat({GroupNames{1}}, numel(g1),1);
g2_idx = repmat({GroupNames{2}},numel(g2),1);
age_idx = cat(1,g1_idx,g2_idx);

try
figure;bar([nanmean(g1);nanmean(g2)]')
catch
end
g1 = g1(:);
g2 = g2(:);

all = cat(1,g1,g2);

if ~ full_levels
lvl_idx = repmat(levels,1,length(all)/length(levels) ) ;
lvl_idx = lvl_idx(:);
else 
    lvl_idx = cat(1,levels(:),levels(:));
end 



[x,main_effects,stats,z]= anovan(all,{age_idx,lvl_idx},'model','interaction','varnames',varnames);

maineffectsplot(all,{age_idx,lvl_idx},'varnames',{'age','level'})
if SaveName 
    suptitle(sprintf('Main Effects for: %s', SaveName'))
    saveas(gcf,sprintf('%s -MainEffects.pdf',SaveName) )
end 

figure
all_stats = multcompare(stats,'Dimension',[1 2]);

if SaveName
    saveas(gcf,sprintf('%s -Comparision.pdf',SaveName) )
end 


