function BD = Band2D(DFF,Lvl)

% takes a df_by_level X neuron object and returns the bandwidth for each
% obect.


m = max(DFF);

DFF = DFF ./ m;
if ~exist('Lvl','var')
    BD = sum(DFF>= .75 );
    return
else
    
   DFF = reshape(DFF,4,8,[]);
    
    for ii =1:Lvl
        BD(ii,:) = squeeze(sum(DFF(ii,:,:) >= .75,2));
    end
end
end