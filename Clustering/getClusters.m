function clusters = getClusters(Passive,expt_idx)

  
idx = any(Passive.experiment_list == find(expt_idx)',2);

clusters = Passive.Clusters.Clusters(idx);