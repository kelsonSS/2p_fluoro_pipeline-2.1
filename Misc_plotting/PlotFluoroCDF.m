function dff = PlotFluoroCDF(DF,type,lvl_idx)

    % type = max or mean
    if ~exist('lvl_idx','var')
         lvl_idx = true(size(DF.DFF,2),1);
    end 

         
    nn_idx = DF.Clean_idx & DF.active{:,2} > 0 ;
    dff = DF.DFF(1:end-1,lvl_idx,nn_idx);
    
    
    
    if strcmp(type,'max')
        dff = squeeze(max(max(dff,[],2))) ;
    else 
        dff = squeeze(nanmean(nanmean(dff,2)));
    end
    
    figure;
    cdfplot(dff)
    title(sprintf('%s Fluorescence',type))
        
