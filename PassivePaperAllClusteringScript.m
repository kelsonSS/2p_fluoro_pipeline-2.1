% FRA-Clustering script 
FRA_noise.K_means_clusters = Cluster_DF(FRA_Noise)
FRA_noise.H_clust = Cluster_DF(FRA_Noise,'H-Clust', 11);
FRA_noise.Old_neurons = ismember(FRA_noise.experiment_list(FRA_noise.Clean_idx),FRA_noise.Old_expt_idx);

FRA_noise.Young_neurons = ~ismember(FRA_noise.experiment_list(FRA_noise.Clean_idx),FRA_noise.Old_expt_idx);

% test to ensure we've captured all the neurons 
size(FRA_noise.K_means_clusters,1)-sum(FRA_noise.Young_neurons) == sum(FRA_noise.Old_neurons); 



% create histograms of which clusters the old and young neurons are part of
PassiveClusterPlot(FRA_noise.K_means_clusters,{FRA_noise.Old_neurons,FRA_noise.Young_neurons})

title('K-means Cluster identity')
legend({'All','Old','Young'},'location','northeast')

f = get(gca(),'Children');
set(f(3),'FaceColor','k')
set(f(2),'FaceColor', [.8,.8,.8] )
set(f(1),'FaceColor', [1,1,1] )
xlabel('Cluster Number')
ylabel('(%)')

PassiveClusterPlot(FRA_noise.H_clusters,{FRA_noise.Old_neurons,FRA_noise.Young_neurons})

title('Heirarchical-Clustering Cluster identity')
legend({'All','Old','Young'},'location','northeast')

f = get(gca(),'Children');
set(f(3),'FaceColor','k')
set(f(2),'FaceColor', [.8,.8,.8] )
set(f(1),'FaceColor', [1,1,1] )
xlabel('Cluster Number')
ylabel('(%)')
%% repeat with FRA_in Quiet
FRA_quiet =Fluoro_to_Table(FRA_list_quiet);
FRA_quiet.K_means_clusters = Cluster_DF(FRA_quiet)
FRA_quiet.H_clusters = Cluster_DF(FRA_quiet,'H-Clust', 12);

FRA_quiet.Old_expt_idx = [10, 11,12,13,14];
FRA_quiet.Old_neurons = ismember(FRA_quiet.experiment_list(FRA_quiet.Clean_idx),FRA_quiet.Old_expt_idx);

FRA_quiet.Young_neurons = ~ismember(FRA_quiet.experiment_list(FRA_quiet.Clean_idx),FRA_quiet.Old_expt_idx);

% test to ensure we've captured all the neurons 
size(FRA_quiet.K_means_clusters,1)-sum(FRA_quiet.Young_neurons) == sum(FRA_quiet.Old_neurons); 



% create histograms of which clusters the old and young neurons are part of
PassiveClusterPlot(FRA_quiet.K_means_clusters,{FRA_quiet.Old_neurons,FRA_quiet.Young_neurons})

title('K-means Cluster identity')
legend({'All','Old','Young'},'location','northeast')

f = get(gca(),'Children');
set(f(3),'FaceColor','k')
set(f(2),'FaceColor', [.8,.8,.8] )
set(f(1),'FaceColor', [1,1,1] )
xlabel('Cluster Number')
ylabel('(%)')

PassiveClusterPlot(FRA_quiet.H_clusters,{FRA_quiet.Old_neurons,FRA_quiet.Young_neurons})

title('Heirarchical-Clustering Cluster identity')
legend({'All','Old','Young'},'location','northeast')

f = get(gca(),'Children');
set(f(3),'FaceColor','k')
set(f(2),'FaceColor', [.8,.8,.8] )
set(f(1),'FaceColor', [1,1,1] )
xlabel('Cluster Number')
ylabel('(%)')
