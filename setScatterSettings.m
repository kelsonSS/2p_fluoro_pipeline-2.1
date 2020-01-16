function setScatterSettings(gcf)

ax = findall(gcf,'Type','Axes');
lgnd = findall(gcf,'Type','Legend');

set(ax,'visible','off')
set(findall(gcf,'Type','text'),'visible','on')
set(ax,'color','none')
set(lgnd,'Location','east')
lgnd_pos = get(lgnd,'position');
% shift right
lgnd_pos = cellfun(@(x) x+[.1 0 0 0] ,lgnd_pos,'UniformOutput',0);
for ii = 1:length(lgnd_pos)
    set(lgnd(ii),'position',lgnd_pos{ii})
end 




