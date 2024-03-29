function PlotCellResponseTiming(DF,lvl)


 % create indicies  
if exist('lvl','var')
    lvl_idx = any(DF.FreqLevelOrder{:,2} == lvl , 2);
    assert(sum(lvl_idx) >0 , 'Incorrect level used') 
else
    lvl_idx = true(1,size(DF.DFF,2));
end 
    
    
  
   active_idx=  DF.active{:,2}>0  ;
     
   % get average reponse
   DFF = squeeze(nanmean(DF.DFF(1:end-2,lvl_idx,active_idx),2));
   
   % get expt_idx 
    expt_list = DF.experiment_list(active_idx);
    
    
    n_expts = max(expt_list);
    
  
    timing_col_id = 1;
    for expt_id = 1:n_expts
        
       DFF_temp = DFF(:, expt_list == expt_id);
   
       if size(DFF_temp,2) < 20 
           continue
       end 
       % get first derivative 
   
   
   DFF_temp = abs(diff(DFF_temp));
   
   [~,timing] =  max(DFF_temp);
   
    timing_prc(timing_col_id,:) = histcounts(timing,[1:9:size(DFF_temp,1)+4])./length(timing);
    timing_col_id = timing_col_id +1;
    end 
    
    try
   [~,~,stats] = anova1(timing_prc, [] ,'off');
  x= multcompare(stats,'Display','off');
  x(1:15,:)  
    catch
    end 
    figure; bar(mean(timing_prc),'BarWidth',1)
    hold on 
    errorbar(mean(timing_prc), std(timing_prc) / sqrt(n_expts) * 1.86 ,'.')
    
           
         