%% Clustering 
% inputs 
clust_num = 10;

% 
DFF2 = squeeze(nanmean(DFF,2);

DFF_ab_max = max(abs(DFF2));
DFF_normalized = DFF2./DFF_ab_max;



Z = linkage(DFF_corr,'Complete','correlation');
clusters = cluster(Z,clust_num);

%%  plotting 
m = max(clusters);
figure
rc = numSubplots(max(clusters));
for clust_iter = 1:m 
    
    % find members of cluster
    DFF_idx = clusters == clust_iter;
    test_DFF = DFF_normalized(:,DFF_idx);
    clust_n = sum(DFF_idx);
    
    % analyze
    mu = nanmean(test_DFF,2);
    sigma = nanstd(test_DFF,[],2);
   
    % plot
    subplot(rc(1),rc(2), clust_iter)
    shadedErrorBar([],mu,sigma)
    title( sprintf('%d Neurons', clust_n) );
end 
    
    