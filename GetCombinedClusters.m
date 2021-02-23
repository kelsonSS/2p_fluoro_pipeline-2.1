function [clust,old_clusters,young_clusters]= GetCombinedClusters(FRA_combined_noise,centroids)

if exist('centroids','var')

    clust = Cluster_DF(FRA_combined_noise,'K-means',20,'normalized',Centroids);
else 
    clust = Cluster_DF(FRA_combined_noise);
    
end 

m = max(clust.Clusters);

old_clust = clust.Clusters(FRA_combined_noise.young_idx ==0);
young_clust = clust.Clusters(FRA_combined_noise.young_idx ==1);

old_clusters = old_clust;
young_clusters = young_clust;

old_clust(old_clust ==0)= [];
young_clust(young_clust ==0)= [];

old_counts =  histcounts(old_clust,m) / length(old_clust);
young_counts =  histcounts(young_clust,m) / length(young_clust);

figure;bar(cat(1,young_counts,old_counts)')
pause




