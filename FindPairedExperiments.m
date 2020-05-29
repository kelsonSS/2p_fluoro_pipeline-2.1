function pairs =  FindPairedExperiments(DF1,DF2)
 
DF1_paths = DF1.DataDirs;
DF2_paths = DF2.DataDirs;

 if isstruct(DF1_paths); DF1_paths =  DF1_paths.DataDirs;end
 
 if isstruct(DF2_paths); DF2_paths =  DF2_paths.DataDirs;end
 
 % find DF pairs
 DF1_paths = cellfun(@fileparts,DF1_paths,'UniformOutput',0);
 DF2_paths = cellfun(@fileparts,DF2_paths,'UniformOutput',0);
 
 path1_idx = [1:numel(DF1_paths)]';
 [bool,path2_idx] = ismember(DF1_paths,DF2_paths);
 
pairs =  cat(2,path1_idx(bool),path2_idx(bool));
 
% clean DFF

clean_idx1 =  ismember(DF1.experiment_list,pairs(:,1));
DF1_Clean = DF1.DFF_norm(:,:,clean_idx1);

clean_idx2 = ismember(DF2.experiment_list,pairs(:,2));
DF2_Clean = DF2.DFF_norm(:,:,clean_idx2);


Cluster1 =  Cluster_DF(DF1_Clean,'H-Clust',6);
Cluster2 =  Cluster_DF(DF2_Clean,'H-Clust',6);



figure
cluster_counts = histcounts2(Cluster2,Cluster1);
normed = cluster_counts./sum(cluster_counts,2);
imagesc(normed);
ylabel('cluster Pre')
xlabel('cluster Post')


 
     