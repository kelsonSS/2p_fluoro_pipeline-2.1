function Plot_ClusterDiversity(Clusters)

% This function takes in a list of integer cluster ids [1,2,1,0,3] and
% returns plots showing the diversity of those clusters 
% note: the program assumes that zero's are entries that were not clustered
% and ignores them 


if isstruct(Clusters)
    Plot_ClusterDiversity(Clusters.Clusters)
end 

m = max(Clusters);

Clusters(Clusters == 0) = [];

counts = histcounts(Clusters,m);
figure
bar(counts / numel(Clusters)) % get Fraction in each cluster

title('Cluster Diversity')
