function plotShadedErrorBar(DFF,color)

% this function is the standard way I call shaded errorbar

% -KS 2021

if ~exist('color','var')
    color = 'k'
end 
if isnumeric(color)
shadedErrorBar([],squeeze(nanmean(DFF,2))...
        ,squeeze(nanstd(DFF,[],2))./ sqrt(size(DFF,2))  * 1.96,{'markerfacecolor',color})
else 
    shadedErrorBar([],squeeze(nanmean(DFF,2))...
        ,squeeze(nanstd(DFF,[],2))./ sqrt(size(DFF,2))  * 1.96,color)
end