function Plot_ClusterDiversityByAnimal(Clusters,expt_ids)


if isstruct(Clusters)
    Plot_ClusterDiversityByAnimal(Clusters.Combined_Classes,Clusters.experiment_list)
    return

end 


m = max(Clusters);

N_expts = max(expt_ids);
Neurons_expt = [];
expt_prc = [];
for expt_idx = 1:N_expts

Clusters_expt = Clusters(expt_ids == expt_idx);

Clusters_expt(Clusters_expt == 0) = [];

neurons = numel(Clusters_expt);
counts = histcounts(Clusters_expt,m);

expt_prc(expt_idx,:) = counts / neurons;
Neurons_expt(expt_idx) = neurons;
end 

expt_prc = expt_prc(Neurons_expt > 20,:);

active = sum(expt_prc > .05,2);
figure 
bar(mean(active));
hold on 
errorbar( mean(active),std(active) ./ sqrt(N_expts) * 1.96 , '.' ) ;
set(gca,'Ylim',[0 m+1])
title('Cluster Diversity')

figure 
bar(mean(expt_prc));
hold on 
errorbar( mean(expt_prc),std(expt_prc) ./ sqrt(N_expts) * 1.96 , '.' ) ;




