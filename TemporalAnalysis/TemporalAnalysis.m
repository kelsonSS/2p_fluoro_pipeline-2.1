function timings = TemporalAnalysis(DF)

    idx = DF.Clean_idx & DF.active{:,2} > 0 ;
    dff = DF.DFF_Z(1:end-1,:,idx);
    
    dff = squeeze(nanmean(dff,2));
 
    
    % find location of maximum value in derivative 
   [m_DFF, m_loc]  = max(diff(dff));
   
   
   timings = histcounts(m_loc,[1:5:size(dff,1)+1]);
   
   timings = timings / size(dff,2);
   
   figure;bar(timings,'BarWidth',1);

   xticklabels(1:1:length(xticks))
   
   
   
   