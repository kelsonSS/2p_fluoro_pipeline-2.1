function Flat= Flatten_Corr(Corr)
% takes the Corr object from Correlations function 
% and flattens it across experiments

Corr = cellfun(@(x) x(:),Corr,'UniformOutput',0);

if size(Corr,2) == 1 ;
    
    Flat = Cell2mat(Corr);
else
    for lvl =  1:size(Corr,2)
    Flat{ii} = cell2mat(Corr(:,lvl));
    end 
end 
    