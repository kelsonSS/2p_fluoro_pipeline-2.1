function dff = PlotFluoroCDF(DF,type,lvl,SaveName)

   
   
   
    % type = max or mean
    if ~exist('lvl','var') || ischar(lvl) || isempty(lvl)
         lvl = true(size(DF.DFF,2),1);
         lvl_str = 'All'
    else 
        lvl_str = num2str(lvl);
 
    lvl = DF.FreqLevelOrder{:,2} == lvl;
    end 
    
    
    if ~exist('SaveName','var')
        SaveName = [];
    end 
    
    
    nn_idx = DF.Clean_idx & DF.active{:,2} > 0 ;
    
 
        
    
    dff = DF.DFF(1:end-1,lvl,nn_idx);
    
    
    
    switch type 
        case 'max'
        dff = squeeze(max(max(dff,[],2))) ;
        case 'mean'
        dff = squeeze(nanmean(max(dff,[],2)));
        case 'min'
        dff = squeeze(min(min(dff,[],2))) ;
        case 'min-mean'
        dff = squeeze(min(min(dff,[],2))) ;

    end
    
    figure;
    cdfplot(dff)
    title_str = sprintf('%sFluorescence-%s',type,lvl_str)
    title(title_str)
    if SaveName
        saveas(gcf, sprintf('%s-%s.pdf', SaveName, title_str) )
    end 