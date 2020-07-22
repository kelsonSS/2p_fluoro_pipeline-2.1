function Plot_Clusters(DF,norm,clusters)
% standalone function to plot clusters from TN object after calling the
% Cluster_DF in the resultant plot cluster object.

if ~exist('norm','var')||isempty(norm) ;norm = 'normalized';end
if ~exist('clusters','var');clusters = DF.Class_idx; end 


switch norm 
    case 'normalized'
        DFF_norm = squeeze(nanmean(DF.DFF_norm,2));
    case 'Z-score'
        DFF_norm = squeeze(nanmean(DF.DFF_Z,2));
    otherwise
        error('norm must be normalized or Z-score')
end 

m = max(clusters);

for clust_num = 1:m
% find members of cluster
    DFF_idx = clusters == clust_num;
    test_DFF = DFF_norm(:,DFF_idx);
   
    
    
    % analyze
    mu = nanmean(test_DFF,2);
    sigma = nanstd(test_DFF,[],2);
    D = nanmean(diff(test_DFF),2);
     clust_n = sum(DFF_idx);
    
     % package
    avg_trace{clust_num}.mu = mu;
    avg_trace{clust_num}.sigma = sigma ;
    avg_trace{clust_num}.clust_n  = clust_n;
    avg_trace{clust_num}.CI = sigma/sqrt(clust_n) * 1.96 ;
    avg_trace{clust_num}.D = D;
    
 
end 


rc = numSubplots(m);
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

    