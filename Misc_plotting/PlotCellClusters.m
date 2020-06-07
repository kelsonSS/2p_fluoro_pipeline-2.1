function PlotCellClusters(Passive,manual_clust_ids,manual_colors)
% This function takes the Passive experiment object with a clusters 
% parameter and outputs one image per experiment showing the locations and
% cluster identies of each cluster 
%
% manual clust_ids - option to use a manual clustering for each experiment
%   %% note- there must be as many ids as there are neurons  

    if exist('manual_colors','var')
        colors = manual_colors;
    else 
        colors = {'k','c','r','b','y'};
    end 

%get cluster list 
    if exist('manual_clust_ids','var')
        cluster_ids = manual_cluster_ids;
        clear manual_clust_ids
    else
        clust_ids = Passive.Class_idx;
    end 

expts = unique(Passive.experiment_list);
for expt_idx = 1:length(expts)
    
   % first we'll load the cell directory and images 
    expt_num =  expts(expt_idx);
    expt_directory  =Passive.DataDirs{expt_num};
    
    
    
    % update progress bar
    if exist('text_length_to_delete','var')
        fprintf(repmat('\b',1,text_length_to_delete));
    end   
    text = sprintf('Running Experiment:%d %d/%d',expt_num,expt_idx,length(expts));
    fprintf(text)
    text_length_to_delete = length(text);
   
    
    % these will load the avg_image as avgGreenFilteredImage
    try
        load(fullfile(expt_directory,'AvgImg.mat'));
   
    catch
        try
        load(fullfile(expt_directory,'AvgGreenImage.mat'));
        catch
        end 
    end 
    
    % then load the cell definition files
      cdef = Passive.CellID{expt_num};
      
     %Find cells belonging to current experiment 
     experiment_cell_idx  = Passive.experiment_list == expt_num;
     classes = Passive.Class_idx(experiment_cell_idx);

    
     
 % plot the average image 
    figure
      title(sprintf('Field #%s', expt_num)); 
     try 
         imagesc(avgGreenfilteredImage)
     catch
         imagesc(zeros(512));% in case no avg image exists 
     end 
    hold on 
    colormap gray
 for current_cluster = 1:max(classes)+1 % classes are 0-indexed 
    
     cells_to_plot = find(classes == current_cluster);
     
 % plot the rings with the appropriate colors 
    for ii = 1:length(cells_to_plot)
        pp = cells_to_plot(ii);
        xpt = cdef.ptsIdx(pp,2);
        ypt = cdef.ptsIdx(pp,3);
        %plot outer edges of rings
        plot(xpt + cdef.smRoiBoundaries{pp}(:,3) .*...
             (cos(cdef.smRoiBoundaries{pp}(:,1))) ,...
        ypt + cdef.smRoiBoundaries{pp}(:,3) .*...
        (sin(cdef.smRoiBoundaries{pp}(:,1))),colors{current_cluster} ) 
        %plot inner edges of rings 
        plot(xpt + cdef.smRoiBoundaries{pp}(:,2) .*...
        (cos(cdef.smRoiBoundaries{pp}(:,1))),...
        ypt + cdef.smRoiBoundaries{pp}(:,2) .*...
        (sin(cdef.smRoiBoundaries{pp}(:,1))),colors{current_cluster} ) 
   
    end  
    
    
    
    
    
   
 end 
end 






  