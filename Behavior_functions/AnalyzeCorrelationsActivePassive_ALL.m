function out = AnalyzeCorrelationsActivePassive_ALL(young,old,fields_to_process)
    
if  ~exist('fields_to_process','var')
    field_names = fieldnames(young);
else 
    field_names = fields_to_process;
end 
    
    
    for fld_idx = 1:length(field_names)
        
        fn = field_names{fld_idx};
        
    out.(fn) = AnalyzeCorrelationsActivePassive(young.(fn).Stats,old.(fn).Stats,fn);
    
    
     plotCorrsActivePassive(young.(fn),old.(fn),fn)
     out.(fn).ActivePassive = getActivePassiveCorrStats(young.(fn),old.(fn),fn);
     
    out.(fn).r_plus = getSICorr(young.(fn),old.(fn),true,[fn '-r_plus']);
    
    out.(fn).r_minus = getSICorr(young.(fn),old.(fn),false,[fn '-r_minus']);
    
    out.(fn).Corr_by_dist = CompareCorrsByDistanceBehavior(young.(fn).Corr_by_distance,...
        old.(fn).Corr_by_distance,['DistanceCorrs-' fn]);
    
    
    %% Hit Vs Early Comparision 
    y_p_20{fld_idx} = young.(fn).SNR20.r_plus';
    y_p_0{fld_idx} = young.(fn).SNR0.r_plus';
    y_m_20{fld_idx} = young.(fn).SNR20.r_minus';
    y_m_0{fld_idx} = young.(fn).SNR0.r_minus';
    

    o_p_20{fld_idx} = old.(fn).SNR20.r_plus';
    o_p_0{fld_idx} = old.(fn).SNR0.r_plus';
    o_m_20{fld_idx} = old.(fn).SNR20.r_minus';
    o_m_0{fld_idx} = old.(fn).SNR0.r_minus';
    
    
    end 
    
    
    
    out.young_r_plus = Compare2AnovaBehavior(y_p_20,y_p_0,...
                          {'young-20','young-0'}, field_names, 'CorrelationsActivePassive_young_r_plus')
    
                   
                      
     out.young_r_minus = Compare2AnovaBehavior(y_m_20,y_m_0,...
                          {'young-20','young-0'}, field_names, 'CorrelationsActivePassive_young_r_plus')
                      
    out.old_r_plus = Compare2AnovaBehavior(o_p_20,o_p_0,...
                          {'old-20','old-0'}, field_names, 'CorrelationsActivePassive_old_r_plus')
    
     out.old_r_minus = Compare2AnovaBehavior(o_m_20,o_m_0,...
                          {'old-20','old-0'}, field_names, 'CorrelationsActivePassive_old_r_minus')
                      
    plotRotated(out.young_r_plus,field_names,'young r+','young_r_plus_rotated')  
    plotRotated(out.young_r_minus,field_names,'young_r-','young_r_minus_rotated')  
    plotRotated(out.old_r_plus,field_names,'old_r+','old_r_plus_rotated')  
    plotRotated(out.old_r_minus,field_names,'old_r-','old_r_minus_rotated')  
  
   
    function plotRotated(data,fields,titleName,SaveName)
        
        
        figure
        PlotGroupedErrorBars(data.means',data.CIs')
        legend(fields)
        xticklabels([20,0])
        xlabel('dB SNR')
        ylabel('Correlation (r)')
        title(titleName)
        
        saveas(gcf,sprintf('%s -ComparisonBars_rotated.pdf',SaveName) )
    
        
        
 
        function out = getSICorr(young,old,above_flag,SaveName)
   
         
    
    
            
            young_data = mungeCorrData(young,above_flag);
            
            
            old_data = mungeCorrData(old,above_flag);
            
            full_save_name = ['CorrelationsActivePassive_' SaveName];
         out = Compare2AnovaBehavior(young_data,old_data,...
                          {'young','old'}, {'20','0'},full_save_name );
        plotRotated(out,{'+20dB','0dB'},SaveName,[full_save_name '_rotated'] )    
                 xticklabels({'Young','Old'})
                 xlabel('Age')
        saveas(gcf, [full_save_name '_rotated.pdf'])
                 
                   
            
            
            
            
            
        
function out = mungeCorrData(d,above_flag)
         
    
    
    if above_flag
       f = 'r_plus'
    else
       f ='r_minus'
    end 
    
    out{1} = d.SNR20.(f)';
    out{2} = d.SNR0.(f)';
    

    function out = getActivePassiveCorrStats(young,old,fn)
               
        
      out.SNR20.r_plus = compareField(young.SNR20,old.SNR20,fn,'20','r_plus_cells');
      out.SNR20.r_minus = compareField(young.SNR20,old.SNR20,fn,'20','r_minus_cells');
      out.SNR0.r_plus = compareField(young.SNR20,old.SNR20,fn,'0','r_plus_cells');
      out.SNR0.r_minus = compareField(young.SNR20,old.SNR20,fn,'0','r_minus_cells');
   
        function out = compareField(young,old,fn,SNR,fieldName)
            Y{1} = young.Passive(young.(fieldName))';
            Y{2} =young.Active(young.(fieldName))';
            O{1} = old.Passive(old.(fieldName))';
            O{2} = old.Active(old.(fieldName))';
            
            SaveName = sprintf('CorrelationsActivePassive_%s_%s_%s',fn,SNR,fieldName);
            out = Compare2AnovaBehavior(Y,O,{'Young','Old'},{'Passive','Active'},SaveName);
        
    function  plotCorrsActivePassive(young,old,SaveName)
       
        
        figure 
        subplot(1,2,1)
        
       plotField(young.SNR20,old.SNR20,'r_plus_cells');
       plotField(young.SNR20,old.SNR20,'r_minus_cells');
        title( '+20dB SNR' )
        
         subplot(1,2,2)
        
         plotField(young.SNR0,old.SNR0,'r_plus_cells');
        plotField(young.SNR0,old.SNR0,'r_minus_cells');
        title( '+0dB SNR' )
        
        saveas(gcf,sprintf('CorrelationsActive_v_Passive-%s-Errorbars.pdf',SaveName))
        
   
        
    function plotField(young,old,fieldname)
        hold on 
        
       mu = [];
       err = [];
      
       
  
           
      %plot young 
       [mu(1),err(1)] = getMeansAndCI(young.Passive(young.(fieldname)) );
       [mu(2),err(2)] = getMeansAndCI(young.Active(young.(fieldname)) );      
       
       errorbar(mu,err,'b','LineWidth',2)
       % plot old
        mu = [];
       err = [];
       [mu(1),err(1)] = getMeansAndCI(old.Passive(old.(fieldname)) );
       [mu(2),err(2)] = getMeansAndCI(old.Active(old.(fieldname)) ); 
        
        errorbar(mu,err,'y','LineWidth',2)
        
        plot([.9 2.1], [0 0],'k--' )
        xlim([.9 2.1])
        
        xticks([1 2])
        ylim([-.15 .4])
        xticklabels({'Passive','Active'})
    
function [mu CI] = getMeansAndCI(data)
    mu = nanmean(data);
    CI = nanstd(data) / sqrt(length(data)) * 1.96;
            
        
    
        
        
        
        
        
   
    