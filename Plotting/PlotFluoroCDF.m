function dff = PlotFluoroCDF(DF,type,lvl)

    % type = max or mean
    if ~exist('lvl','var')
         lvl = true(size(DF.DFF,2),1);
    end 
    
    lvl = DF.FreqLevelOrder{:,2} == lvl;
    
    assert( sum(lvl) > 0)
         
    nn_idx = DF.Clean_idx & DF.active{:,2} > 0 ;
    dff = DF.DFF(1:end-1,lvl,nn_idx);
    
    
    
    if strcmp(type,'max')
        dff = squeeze(max(max(dff,[],2))) ;
    else 
        dff = squeeze(nanmean(nanmean(dff,2)));
    end
    
    figure;
    cdfplot(dff)
    title(sprintf('%s Fluorescence',type))