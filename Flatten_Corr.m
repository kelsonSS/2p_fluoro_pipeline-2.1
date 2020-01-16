function Flat= Flatten_Corr(Corr)
% takes the Corr object from Correlations function 
% and flattens it across experiments

Corr = cellfun(@(x) x(:),Corr,'UniformOutput',0);

if isrow(Corr)|| iscolumn(Corr) ;
    try
    Flat = cell2mat(Corr);
    catch
    Flat = cell2mat(Corr');
    end 
else
    for lvl =  1:size(Corr,2)
    Flat{lvl} = cell2mat(Corr(:,lvl));
    end 
end 
    