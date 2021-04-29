function [all_stats,main_effects] = CompareOldYoung(old,young,levels,full_levels)

% this function will create an anova to compare two groups
% at different levels 
%
% input 
% old/young  matrices in a samples X levels format
% levels - 1 x N cell of different levels used {'90';'20';'10;'0'}

if ~exist('full_levels','var')
    full_levels = false;
end 


old_idx = repmat({'old'}, numel(old),1);
young_idx = repmat({'young'},numel(young),1);
age_idx = cat(1,old_idx,young_idx);

try
figure;bar([nanmean(old);nanmean(young)]')
catch
end
old = old(:);
young = young(:);

all = cat(1,old,young);

if ~ full_levels
lvl_idx = repmat(levels,1,length(all)/length(levels) ) ;
lvl_idx = lvl_idx(:);
else 
    lvl_idx = cat(1,levels(:),levels(:));
end 



[x,main_effects,stats,z]= anovan(all,{age_idx,lvl_idx},'model','interaction');

all_stats = multcompare(stats,'Dimension',[1 2]);

