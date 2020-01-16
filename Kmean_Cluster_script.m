%%Cluster Script
% this script uses an optimized version of Kmeans clustering 
% to find the ideal number of clusters and then display the average 
% FRA plot of each group

% This plot will first diplay a number of Kmeans plots determined by k
% below. The final plot will be the plot determined by the elbow method to
% be the optimal number of groups

%% Cluster creation 
for k = 5:20
    cluster_km = kmeans(df_by_level,k);% df_norm or df_by_level
    
    m = max(cluster_km);
    df_clstr = zeros(m,size(df_by_level,2));
    df_clstr_n = zeros(m,1);
    
    for grp =  1:m
        % extract group id
        grp_idx = (cluster_km == grp);
        grp_n =  sum(grp_idx);
        
        
        % create average DFF plot
        g_tmp =  df_by_level(grp_idx,:);
        g_tmp = nanmean(g_tmp);
        
        df_clstr(grp,:) = g_tmp;
        df_clstr_n(grp,:) = grp_n;
        
        
    end
    
    
    %% Cluster Plotting
    m  = max(cluster_km);
    rc = numSubplots(size(df_clstr,1));
    figure
    for ii = 1:m
        
        subplot(rc(1),rc(2),ii)
        axis off
        
        myFRA(uFreqs,uLevels',df_norm(ii,:),ii,false)
        title(sprintf('Group %d: %d Neurons',ii,df_clstr_n(ii)))
        
        
    end
end


%Optimal plotting using k_means_opt
cluster_km = kmeans_opt(df_by_level);% df_norm or df_by_level

m = max(cluster_km);
df_clstr = zeros(m,size(df_by_level,2));
df_clstr_n = zeros(m,1);

for grp =  1:m
    % extract group id 
    grp_idx = (cluster_km == grp);
    grp_n =  sum(grp_idx);
    
    
    % create average DFF plot
    g_tmp =  df_by_level(grp_idx,:);
    g_tmp = nanmean(g_tmp);
    
    df_clstr(grp,:) = g_tmp;
    df_clstr_n(grp,:) = grp_n; 




end
rc = numSubplots(size(df_clstr,1));
figure
for ii = 1:m  
   
   subplot(rc(1),rc(2),ii)
   axis off 
   
   myFRA(uFreqs,uLevels',df_norm(ii,:),ii,false)
   title(sprintf('Group %d: %d Neurons',ii,df_clstr_n(ii)))
   
    
end