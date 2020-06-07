function [out_clusters,neuron_id,avg_trace,...
          var_explained,cluster_centroids] = Cluster_DF(DF,varargin) 
%%% [clusters,DF] = Cluster_DF(DF,varargin) returns the cluster identity
%%% [clusters] and average cluster Fluorescence (avg_trace) for each cluster.
%%% Specifically, avg_trace is a cell containing the mean and STD of each cluster
%%% This will also create a plot of the average fluorescence traces. 
%%% Clustering will be performed either using kmeans_opt function of by
%%% heiraricha clustering, depening on which algo is selected. 
%%% The input DF is either:
%%%      a time X trial X neuron  matrix 
%%% or 
%%%      the DF output of Fluoro_to_Table 
%%%                                                                       
%%%   Cluster_DF(DF) returns the optimal number of clusters via KMEANS_OPT
%%%
%%%   Cluster_DF(DF,MODE) will allow the user to use heirarical clustering
%%%   or kmeans (arguments, 'K-means', 'H-Clust' default ='K-means')
%%%   Cluster_DF(DF,MODE,MAX_CLUST) will allow the user to define the
%%%   number of clusters (default=20)
%%%   
%%%   Kelson Shilling-scrivo 2020 

%% Clustering 
% inputs 

if nargin>1, mode = varargin{1}; else, mode = 'K-means';end
if nargin>2, max_clust = varargin{2}; else, max_clust = 20;end
if nargin>3, norm_mode = varargin{3};else, norm_mode = 'normalized';end    
if nargin>4, centroids = varargin{4};end ;

if isstruct(DF)
    DF.DFF2 = squeeze(nanmean(DF.DFF,2));
else
    temp = DF;
    clear DF
    DF = struct('DFF',temp, 'DFF2', squeeze(nanmean(temp,2)) )
    clear temp
end

 %%  normalizing to absolute max 
switch norm_mode
    case 'normalized'
        DFF_ab_max = max(abs(DF.DFF2));
         DF.DFF_norm = DF.DFF2./DFF_ab_max;

%% normalizing to baseline corrected z-score (decide which one)
    case 'Z-score'
        DF.DFF_norm = squeeze(nanmean(DF.DFF_Z,2));
      
    otherwise
        error('norm_mode must be set to normalized or Z-score')
end 

  neuron_id =1:size(DF.DFF_norm,2);

% prealloate output structure
out_clusters = [1:length(neuron_id)]';
out_clusters(:,2) = 0;


% clean and  subscript neurons 
if isfield(DF,'Clean_idx') & isfield(DF,'active')
     active_idx = DF.active{:,2} > 0;
     neuron_id = 1:size(DF.DFF_norm,2);
        neuron_id = neuron_id(active_idx & DF.Clean_idx);
     DF.DFF_norm = DF.DFF_norm(:,active_idx & DF.Clean_idx );
     
  
end

% Create linkage plot and clusters by looking at the correlations between
% neurons

 var_explained = [];
switch mode
   
    case 'K-means'
        [clusters,cluster_centroids,...
         ~,~,var_explained] = kmeans_opt(DF.DFF_norm',max_clust);
    case 'Paired-K-means'
        [clusters,cluster_centroids] = kmeans(DF.DFF_norm',max_clust,'Start',centroids);
     
    case 'H-Clust'
        DFF_corr = corr(DF.DFF_norm);    
        Z = linkage(DFF_corr,'Complete','correlation');
        clusters = cluster(Z,'MaxClust',max_clust);
    otherwise
        error('no such mode exists!, use K-means or H-Clust')
end

%% Sorting 

% presort clusters by waveform activity 
try
    m = max(clusters(:));
    rc = numSubplots(m); % rc is needed for plotting later 
catch 
    m = max_clust;
    rc = numSubplots(max_clust); % see numsubplot for more info
end 
% find waveform dymanics and save them in avg_trace cell
for clust_num = 1:m
% find members of cluster
    DFF_idx = clusters == clust_num;
    test_DFF = DF.DFF_norm(:,DFF_idx);
    avg_trace{clust_num}.clust_n = sum(DFF_idx);
    
    % analyze
    mu = nanmean(test_DFF,2);
    sigma = nanstd(test_DFF,[],2);
    D = nanmean(diff(test_DFF),2);
   
    avg_trace{clust_num}.mu = mu;
    avg_trace{clust_num}.sigma = sigma;
    avg_trace{clust_num}.D = D;
 
end 
%% split into positive and negative

% find time in which max value occured for each clluster 
max_times = cellfun(@(x) find(abs(x.mu) == max(abs(x.mu))),avg_trace,'UniformOutput',1);

% find if that value was positive or negative at max
for clust_num = 1:m
   pos_neg(clust_num) =  sign(avg_trace{clust_num}.mu(max_times(clust_num)));
end 

% add sign and sort in descending order (positive first)

[~, clust_order] =  sort(max_times .*pos_neg,'descend');
%% replace cluster numbers to sort everything in new order 
[~,transition_mat] = sort(clust_order);

[C,~,clust_vals] =unique(clusters);
 new_clusters =  transition_mat(clust_vals);
 new_clusters = reshape(new_clusters,size(clusters));

 % populate outputs with new clusters 
out_clusters(neuron_id,2) = new_clusters;
% for compatibility we only use the clusters and not the neuron ID
out_clusters = out_clusters(:,2); 


% sort avg traces to correct order 
avg_trace = avg_trace(clust_order);

clear clust_vals
clear C 




%%  plotting
figure
for clust_num = 1:m 
    

    subplot(rc(1),rc(2), clust_num)
    shadedErrorBar([],avg_trace{clust_num}.mu,avg_trace{clust_num}.sigma);
    title( sprintf('%d Neurons', avg_trace{clust_num}.clust_n) );
    %ylim([-1 1])
     axis tight 
    aa = axis;
    hold on
    
    
    if size(test_DFF,1) > 90
        style = 'noise';
    else 
        style = 'tones';
    end 
    
    style_cluster_plot(aa,style)
end 
    
figure

for clust_num = 1:m 
    

    subplot(rc(1),rc(2), clust_num)
    plot(avg_trace{clust_num}.D,'k','LineWidth',1.5);
    title( sprintf('%d Neurons', avg_trace{clust_num}.clust_n) );
    aa = axis;
    hold on
    
    
    if size(test_DFF,1) > 90
        style = 'noise';
    else 
        style = 'tones';
    end 
    
    style_cluster_plot(aa,style)
end 
    






function style_cluster_plot(aa,style) 
        
    
    plot([aa(1) aa(2)], [0 0 ] ,'k--')
    if style == 'noise'
        plot([30 30], [aa(3) aa(4)],'r--')
         plot([60 60], [aa(3) aa(4)],'g--')
        plot([90 90], [aa(3) aa(4)],'g--')
        plot([120 120], [aa(3) aa(4)],'r--')
        axis tight
        xlim([0 150])
    else 
        plot([30 30], [aa(3) aa(4)],'g--')
        plot([60 60], [aa(3) aa(4)],'g--')
        axis tight
        xlim([0 90])
    end 

    

    