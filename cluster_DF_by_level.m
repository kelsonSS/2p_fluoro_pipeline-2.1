function out = Cluster_DF_by_level(DF,max_clust,clust_style)
% this function extends Cluster_DF by extending it by level
% we will break down each level neuron combination and then loop over to
% find how each neuron changes groups across levels

if ~exist('max_clust');max_clust = 20;end 

if ~exist('clust_style','var');clust_style = 'pooled';end 


handles = DF.handles{1};
tempDF.handles = handles;
tempDF.active = DF.active;


% clean DFF to contain only clean and active cells 
clean_idx = DF.Clean_idx & ( DF.active{:,2}>0 );
DF.DFF = DF.DFF(:,:,clean_idx);
DF.DFF_Z = DF.DFF_Z(:,:,clean_idx);
cleanNeuronNumber = sum(clean_idx);

%% replace any Nans with 0s - testing
DF.DFF(isnan(DF.DFF)) = 0;

if isfield(DF,'DFF_Z')
    DF.DFF_Z(isnan(DF.DFF_Z)) = 0;
end 


%% main loop for extracting levels over experiments

%todo-  all tones and levels are in same positions they should be
%so this loop could be reduced to a selection but would be less general
tempDF.DFF = [];
tempDF.DFF_Z =[];
Levels= DF.FreqLevelOrder{:,2};
uLevels = unique(Levels);
for lvl = 1:length(uLevels)
    out.levels(lvl) = uLevels(lvl);
    lvl_DF.DFF =[];
    lvl_DF.DFF_Z =[];
    for expts = 1:length(handles)
        handles = DF.handles{expts};
        lvl_DF.DFF = cat(3, lvl_DF.DFF,...
                     DF.DFF(:,Levels == uLevels(lvl),: ));
        lvl_DF.DFF_Z = cat(3, lvl_DF.DFF_Z,...
                     DF.DFF_Z(:,Levels == uLevels(lvl),: ));
    end
    
    switch  clust_style 
        case 'seperate' % each level has unique clusters 
       [ out.Clusters{lvl}, out.DF{lvl} ]= Cluster_DF(lvl_DF,'K-means',max_clust);
       figure
       title(sprintf('Clusters for level %d', handles.uLevels(lvl)))
        case 'pooled' % finds global clusters and cluster
                      % transisitons across levels

         tempDF.DFF = cat(3,tempDF.DFF,lvl_DF.DFF);
         tempDF.DFF_Z = cat(3,tempDF.DFF_Z,lvl_DF.DFF_Z);
        otherwise 
            error('clust_style must either = pooled or seperate')
    end
    
 
  
  
    
end 
 
% pooled analysis
if strmatch(clust_style,'pooled')
    [out.Clusters, out.DF] = Cluster_DF(tempDF,'K-means',max_clust);

        %% pooled
        % permute into level X neuron order
        out.Clusters = reshape(out.Clusters,length(handles.uLevels),cleanNeuronNumber);
        plot_clusters = out.Clusters;
elseif strmatch(clust_style,'seperate')
    
    
        plot_clusters = cell2mat(out.Clusters)';
       
end 
      
        
        
        num_clust = max(plot_clusters(:));
        % plotting
        PlotPooledClusters(plot_clusters,num_clust,handles.uLevels)

        
 

%% repeat for heirarchical clustering
%[out.H_Clusters,out.DF] = Cluster_DF(tempDF,'H-Clust',num_clust);
%out.H_Clusters = reshape(out.H_Clusters,length(handles.uLevels),cleanNeuronNumber);
%Plotting

%PlotClusters(out.H_Clusters,num_clust,handles.uLevels)


function PlotPooledClusters(C,num_clust,uLevels)
colormap parula
 
 if size(C,1)< 2 ;return; end ;
 
 % convert notations anything above 90db is inf snr
uLevels(uLevels>90) = inf;
 

%

for lvl = size(C,1):-1:2
    current_level = uLevels(lvl);
    next_level = uLevels(lvl-1);
    % create transition matrix between two levels
    
   
    % count transition frequencies
    trans_mat = histcounts2(C(lvl-1,:),C(lvl,:))';
        %% normalize matrix 
        % each column should add to probability 1 
        trans_mat = trans_mat./sum(trans_mat,2);
        % plotting 
        figure
        imagesc(trans_mat)
        % beautification 
        xlabel('End Cluster #')
        ylabel('Starting Cluster #')
        colorbar
        title(sprintf('Transition matrix between %d db and %d db',uLevels(lvl),uLevels(lvl-1) ));
        
        
        
        
        
        
    
    
    
    end 
    
    
 








