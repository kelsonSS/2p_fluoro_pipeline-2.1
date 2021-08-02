function [Out] = CompareOldYoung(old,young,GroupNames,Levels,SaveName,full_levels)

% 
 if ~exist('SaveName','var')
     SaveName = [];
 end 

s1_animal = old{end}
s2_animal = young{end}
 
[Out.Cell.all_stats,Out.Cell.main_effects] = Compare2Anova(old{1},young{1},GroupNames,...
           Levels,[SaveName '-Cell']);
       
[Out.Animal.all_stats,Out.Animal.main_effects] = Compare2Anova(old{end},young{end},GroupNames,...
           Levels,[SaveName '-animal']);       


title_str = sprintf('%s-byAnimal-scatter',SaveName)
       
 PlotScatterLevels(GroupNames,Levels,title_str,s1_animal,s2_animal)
% this function will create an anova to compare two groups
% at different levels 
% 
% input 
% old/young  matrices in a samples X levels format
% levels - 1 x N cell of different levels used {'90';'20';'10;'0'}
% 
% if ~exist('full_levels','var')
%     full_levels = false;
% end 
% 
% 
% old_idx = repmat({'old'}, numel(old),1);
% young_idx = repmat({'young'},numel(young),1);
% age_idx = cat(1,old_idx,young_idx);

% try
% figure;bar([nanmean(old);nanmean(young)]')
% catch
% end
% old = old(:);
% young = young(:);
% 
% all = cat(1,old,young);
% 
% if ~ full_levels
% lvl_idx = repmat(levels,1,length(all)/length(levels) ) ;
% lvl_idx = lvl_idx(:);
% else 
%     lvl_idx = cat(1,levels(:),levels(:));
% end 
% 
% 
% 
% [x,main_effects,stats,z]= anovan(all,{age_idx,lvl_idx},'model','interaction','varnames',{'age','level'});
% 
% maineffectsplot(all,{age_idx,lvl_idx},'varnames',{'age','level'})
% if exist('SaveName','var')
%     suptitle(sprintf('Main Effects for: %s', SaveName'))
%     saveas(gcf,sprintf('%s -MainEffects.pdf',SaveName) )
% end 
% 
% figure
% all_stats = multcompare(stats,'Dimension',[1 2]);
% 
% if exist('SaveName','var')
%     saveas(gcf,sprintf('%s -Comparision.pdf',SaveName) )
% end 
% 

