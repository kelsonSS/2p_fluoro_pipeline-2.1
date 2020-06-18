function out = PlotCellClusters(Passive,manual_clust_ids,manual_colors)
% This function takes the Passive experiment object with a clusters 
% parameter and outputs one image per experiment showing the locations and
% cluster identies of each cluster 
%
% manual clust_ids - option to use a manual clustering for each experiment
%   %% note- there must be as many ids as there are neurons  


%TODO- Currently an issue where the rings are not stored correctly from the
%CellIDS celldefinitionGUI plotting with just cell centers for now. would
%need to regenerate missing roi boundaries to fix this issue. 

    if exist('manual_colors','var')
        colors = manual_colors;
    else 
        % colors = {'k','r','o','y','c','b'};for reference 
        colors = {[1 1 1 ],[1 0 0],[1 .5 0], [1 1 0],[0 1 1] ,[0 0 1]};
   
        
    end 

%get cluster list 
    if exist('manual_clust_ids','var')
        cluster_ids = manual_cluster_ids;
        clear manual_clust_ids
    else
        clust_ids = Passive.Class_idx;
    end 
m = max(clust_ids);
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
   
      % then load the cell definition files
      cdef = Passive.CellID{expt_num};
      
     %Find cells belonging to current experiment 
     experiment_cell_idx  = Passive.experiment_list == expt_num;
     classes = Passive.Class_idx(experiment_cell_idx);

    
     xpts_all = cdef.ptsIdx(:,2);
     ypts_all = cdef.ptsIdx(:,3);
     
     if length(xpts_all)> length(classes)
         xpts_all = xpts_all(1:length(classes));
         ypts_all = ypts_all(1:length(classes));
     elseif length(xpts_all) < length(classes)
         classes = classes(1:length(xpts_all));
     end 
  %KNN Analysis
  responsive_idx = classes ~=0;
  classes = classes(responsive_idx);
  xpts_all = xpts_all(responsive_idx);
  ypts_all = ypts_all(responsive_idx);
 [out.KNNs{expt_num},out.Distances{expt_num} ] = knnsearch(xpts_all,...
                                               ypts_all,'K',10);
  out.KNN_classes{expt_num} =  classes(out.KNNs{expt_num});
   
       
    
    
    % these will load the avg_image as avgGreenFilteredImage
    try
        load(fullfile(expt_directory,'AvgImg.mat'));
   
    catch
        try
        load(fullfile(expt_directory,'AvgGreenImage.mat'));
        catch
            % create an average from raw image 
%             f = fullfile(expt_directory,'GreenChannelRegistered.mat')
%             
        end 
    end 
    
  
   
    
     
     
 % plot the average image 
    figure
      title(sprintf('Field #%s', expt_num)); 
     try 
         imagesc(avgGreenfilteredImage)
         clear avgGreenfilteredImage
     catch
         imagesc(zeros(512));% in case no avg image exists 
     end 
    hold on 
    colormap gray
          title_text = sprintf('Experiment %d, %s',...
              expt_idx,Passive.Expt_names{expt_num});
    title(title_text,'Interpreter','None');
    
 for current_cluster = 1:max(classes) % classes are 0-indexed 
    
     cells_to_plot = find(classes == current_cluster);
     
     cells_to_plot(cells_to_plot>length(cdef.ptsIdx)) =[];
     
     xpts = cdef.ptsIdx(cells_to_plot,2);
     ypts = cdef.ptsIdx(cells_to_plot,3);
     
     scatter(xpts,ypts,'MarkerEdgeColor', colors{current_cluster});
     
     
 % plot the rings with the appropriate colors - currently disabled due to
 %                                              issue with cdef file having 
 %                                               missing smBoundaries 
%     for ii = 1:length(cells_to_plot)
%         pp = cells_to_plot(ii);
%         xpt = cdef.ptsIdx(pp,2);
%         ypt = cdef.ptsIdx(pp,3);
%         %plot outer edges of rings
%         plot(xpt + cdef.smRoiBoundaries{pp}(:,3) .*...
%              (cos(cdef.smRoiBoundaries{pp}(:,1))) ,...
%         ypt + cdef.smRoiBoundaries{pp}(:,3) .*...
%         (sin(cdef.smRoiBoundaries{pp}(:,1))),colors{current_cluster} ) 
%         %plot inner edges of rings 
%         plot(xpt + cdef.smRoiBoundaries{pp}(:,2) .*...
%         (cos(cdef.smRoiBoundaries{pp}(:,1))),...
%         ypt + cdef.smRoiBoundaries{pp}(:,2) .*...
%         (sin(cdef.smRoiBoundaries{pp}(:,1))),colors{current_cluster} ) 
%    
%     end  
%     
%     
    
% gather KNN_classes 

Histc= histcounts(...
    out.KNN_classes{expt_num}(cells_to_plot,:),...
    [.5:1:m+.5]) % edges 
% normalize 
 out.Histc{expt_num,current_cluster}= Histc ./ sum(Histc);
 
 % edge case for no neurons in cluster 
if  isempty(out.Histc{expt_num,current_cluster})||...
        any(isnan(out.Histc{expt_num,current_cluster}))
    out.Histc{expt_num,current_cluster}= zeros(1,m);
end 
        

 end 
 
% 
%  figure
%   bar(cell2mat(out.Histc(expt_num,:)'),'stacked')
%  
 end



Sum_Histc = cell(1,m);

for expt_idx = 1:length(out.Histc)
    for class_idx = 1:m;
        if size(out.Histc{1,class_idx},2) == m
            
            if isempty(Sum_Histc{1,class_idx})
               Sum_Histc{1,class_idx} = out.Histc{1,class_idx}
            else 
                Sum_Histc{1,class_idx} = ...
                    Sum_Histc{1,class_idx} +...
                    out.Histc{1,class_idx};
            end 
        end 
    end 
end 
% for ii = 1: m            
%     figure
%         bar(Sum_Histc{ii})
%     title(sprintf('KNN of cluster %d',ii));
% end 
%    
end 






  