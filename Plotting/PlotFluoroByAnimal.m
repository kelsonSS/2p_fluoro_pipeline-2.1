function results_med = PlotFluoroByAnimal(DF,type,lvl,young_flg)
% type = max or mean
   if ~exist('young_flg','var')
       young_flg = 1;
   end 
   
   if young_flg
       plt_color = 'b';
   else
       plt_color ='k';
   end 
   


    if ~exist('lvl','var')
         lvl_idx = true(size(DF.DFF,2),1);
    else 
    
    lvl_idx = DF.FreqLevelOrder{:,2} == lvl;
    
    assert( sum(lvl_idx) > 0)
    end 
         
    nn_idx = DF.Clean_idx & DF.active{:,2} > 0 ;
    
    n_expts = max(DF.experiment_list);
    
    
   
    

results = {};
%figure
%hold on

expt_ids_plotting ={};
for expt_id = 1:n_expts
       
       
       DFF_temp = DF.DFF(:,lvl_idx, nn_idx & (DF.experiment_list == expt_id) );
       
       
         if strcmp(type,'max')
        results{expt_id} = squeeze(max(max(DFF_temp,[],2))) ;
    else 
        results{expt_id} = squeeze(nanmean(max(DFF_temp,[],2)));
         end
         
    expt_ids_plotting{expt_id} = repmat(expt_id, length( results{expt_id} ), 1 ); 
         
   % cdfplot(results{expt_id})
end


%figure
%boxplot(cell2mat(results'),cell2mat(expt_ids_plotting'))
%title(sprintf( '%s Fluorescence', type))

figure

results_med = cellfun(@median, results );
%boxplot( results_med)
scatter( ones(length(results_med),1) - young_flg , results_med,['.', plt_color])
hold on
axis tight
ylim([ 0, max(ylim)] )
xlim([-1 , 2])




