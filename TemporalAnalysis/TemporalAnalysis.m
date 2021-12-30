function Timings = TemporalAnalysis(DF,lvl,SaveName)


if ~exist('SaveName','var')
    SaveName = ''
end 


bin_size = 10;



 % create indicies  
if exist('lvl','var') && ~strcmp('all',lvl)
    lvl_idx = any(DF.FreqLevelOrder{:,2} == lvl , 2);
    assert(sum(lvl_idx) >0 , 'Incorrect level used') 
else
    lvl_idx = true(1,size(DF.DFF,2));
end 
    

 
  
   active_idx=  DF.active{:,2}>0  ;
     
   % get average reponse
   DFF = squeeze(nanmean(DF.DFF(1:end-2,lvl_idx,active_idx),2));
    
    % get bins 
   hist_bins = [1:bin_size:size(DFF,1)+4];


   % get expt_idx 
    expt_list = DF.experiment_list(active_idx);
    
    

       % get first derivative 
   
   
   DFF = diff(DFF);
   
   [dff_max,timing] =  max(DFF);
   



    exc_idx = dff_max > 0;
    timing = timing(exc_idx);
   

    timing_prc = histcounts(timing,hist_bins)./length(timing);

   



Timings.Cell_Proportion = timing_prc;

    %% plotting
    
    figure; bar(timing_prc,'BarWidth',1)
    hold on 
%   errorbar(mean(timing_prc), std(timing_prc) / sqrt(n_expts) * 1.96 ,'.')
    
    % plotting individual animals with scatter plot
    
    %make a matrix denoting each timebin 
%     
%     n_timebins = size(timing_prc,2);
%     n_expts = size(timing_prc,1);
%     
%     timing = repmat( [1:n_timebins],n_expts,1);
%     
%  scatter(timing(:),timing_prc(:),15, 'MarkerEdgeColor','k','MarkerFaceColor','w')
 ylim([0 .5])
 ticks = yticks;
 yticklabels(ticks * 100)
 ylabel('%')
title([SaveName '-Timing-Cell'],'Interpreter','None')
 if SaveName
     saveas(gcf,sprintf('%s-TimingByCell.pdf', SaveName))
 end
    

           
         
