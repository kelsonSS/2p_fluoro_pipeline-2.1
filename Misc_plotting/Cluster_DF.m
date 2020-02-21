function [clusters,DF] = Cluster_DF(DF,clust_num) 

%% Clustering 
% inputs 
if ~exist('clust_num','var')
    clust_num = 5;
end 
% 
if isstruct(DF)
DF.DFF2 = squeeze(nanmean(DF.DFF,2));
else 
temp = DF;
clear DF
DF = struct('DFF',temp, 'DFF2', squeeze(nanmean(temp,2)) )
clear temp
end

DFF_ab_max = max(abs(DF.DFF2));
DF.DFF_normalized = DF.DFF2./DFF_ab_max;
DFF_corr = corr(DF.DFF_normalized);


Z = linkage(DFF_corr,'Complete','correlation');
clusters = cluster(Z,clust_num);

%%  plotting 
m = max(clusters);
figure
rc = numSubplots(max(clusters));
for clust_iter = 1:m 
    
    % find members of cluster
    DFF_idx = clusters == clust_iter;
    test_DFF = DF.DFF_normalized(:,DFF_idx);
    clust_n = sum(DFF_idx);
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    % analyze
    mu = nanmean(test_DFF,2);
    sigma = nanstd(test_DFF,[],2);
   
    % plot
    subplot(rc(1),rc(2), clust_iter)
    shadedErrorBar([],mu,sigma)
    title( sprintf('%d Neurons', clust_n) );
    ylim([-1 1])
    aa = axis;
    hold on
    
    
    plot([aa(1) aa(2)], [0 0 ] ,'k--')
    plot([30 30], [aa(3) aa(4)],'r--')
    plot([60 60], [aa(3) aa(4)],'g--')
    plot([90 90], [aa(3) aa(4)],'g--')
    plot([120 120], [aa(3) aa(4)],'r--')
    axis tight 
    xlim([0 150])
end 
    
    