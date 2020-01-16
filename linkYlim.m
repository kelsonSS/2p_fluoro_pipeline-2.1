function linkYlim(gcf)
 % this function takes a figure handle with multiple subplots and links the
 % graphs such that they have the same Y-limits as defined by the largest
 % and smallest Y-limts found in the individual subplots

 ax =  findall(gcf,'Type','Axes');

 yl = get(ax,'ylim');
 
 ymin = min(cellfun(@min, yl)) ;
 ymin = ymin+ymin/10;
 
 ymax = max(cellfun(@max, yl)) ;
 ymax = ymax- ymax/10;
 
 set(ax,'ylim',[ymin ymax])
 
 
