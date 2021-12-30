function [Out,max_dist] = getCellDistances(DF)
% this function produces a matrix of all cell distances between neurons 

Out = {};
max_dist  = 0;
cell_flag = 0;

if isstruct(DF)
    n_expts = length(DF.CellID);
    CellIDs = DF.CellID;
    active = DF.active{:,2};
    Expt_list = DF.experiment_list;
elseif iscell(DF)
    n_expts = length(DF)
    cell_flag = 1
else
error('entry should be struct or cell');
end 



for expt_idx =1:n_expts

if cell_flag
    CellID_expt = DF{expt_idx}.CellID{1};
    active = DF{expt_idx}.active{:,2};
    Expt_list = DF{expt_idx}.experiment_list * expt_idx ;

else 
    CellID_expt = CellIDs{expt_idx}
end

microns_px = CellID_expt.micronsPerPixelX;
cell_id_microns = CellID_expt.ptsIdx(:,2:3) * microns_px;

active_idx = active(Expt_list == expt_idx)> 0 ;

cell_id_microns = cell_id_microns(active_idx,:);

dist = squareform(pdist(cell_id_microns));

max_dist_animal = nanmax(dist(:));

max_dist = max( [max_dist,max_dist_animal]);
 
dist(dist == 0) = nan;

Out{expt_idx} = dist;

end