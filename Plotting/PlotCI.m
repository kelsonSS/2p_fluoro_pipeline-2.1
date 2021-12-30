function PlotCI(m,fig_handle,labels)


if ~exist('fig_handle','var')
figure
else
figure(fig_handle)    
end
hold on
bar(nanmean(m))
errorbar([],nanmean(m) ,nanstd(m) / sqrt(size(m,1)) * 1.96,'k.')

if exist('labels','var')
    xticklabels(labels)
end 
    

