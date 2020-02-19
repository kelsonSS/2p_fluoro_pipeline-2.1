function plotAllDFFs(DFF)
    
if isstruct(DFF)
    clean_idx = DFF.Clean_idx;
    DFF = DFF.DFF(:,:,clean_idx);
   

end 

figure
hold on 

for nn = 1:size(DFF,3)
    shadedErrorBar([],squeeze(nanmean(DFF(:,:,nn),2))...
        ,squeeze(nanstd(DFF(:,:,nn),[],2))./320 * 1.96 )
    
    
    
    
end 