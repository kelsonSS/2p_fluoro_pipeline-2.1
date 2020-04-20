function [CLUST,edges] = PassiveClusterPlot(clusters,subsets,edges)
%%% this function takes in the cluster output of cluster-DF which is an
%%% integer list of clusters from 1:N_clusters and will produce a graph
%%% showing the relative population of each cluster


 [CLUST,edges] =  histcounts(clusters,'BinMethod','integer');
 
 CLUST = CLUST./length(clusters);
 
 if exist('subsets','var')
     for s_idx = 1:size(subsets,2)
         
         sub_clust = clusters(subsets{s_idx});
         temp_clust = histcounts(sub_clust ,'BinEdges', edges);
                              
         temp_clust = temp_clust ./ length(sub_clust);
         
         CLUST = cat(1,CLUST,temp_clust)
         
     end 
 end 
 
 figure()
 bar(CLUST')
 
         

 
 
 
 

 
