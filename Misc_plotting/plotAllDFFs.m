function plotAllDFFs(DFF,norm)
    
if isstruct(DFF)
    if isfield(DFF,'Clean_idx')
        clean_idx = DFF.Clean_idx;
    elseif isfield(DFF,'ArtifactIndex')
        clean_idx = ~DFF.ArtifactIndex;
    end 
    
    
    DFF = DFF.DFF(:,:,clean_idx);
end 

figure
hold on 

for nn = 1:size(DFF,3)
    shadedErrorBar([],squeeze(nanmean(DFF(:,:,nn),2))...
        ,squeeze(nanstd(DFF(:,:,nn),[],2))./320 * 1.96 )
    
    
    
    
end 