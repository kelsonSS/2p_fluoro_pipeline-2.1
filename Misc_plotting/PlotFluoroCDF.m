function PlotFluoroCDF(DF,type)

    % type = max or mean

    idx = DF.Clean_idx & DF.active{:,2} > 0 ;
    dff = DF.DFF(1:end-1,:,idx);
    
    if strcmp(type,'max')
        dff = squeeze(max(nanmean(dff,2))) ;
    else 
        dff = squeeze(nanmean(nanmean(dff,2)));
    end
    
    figure;
    cdfplot(dff)
    title(sprintf('%s Fluorescence',type))
        
