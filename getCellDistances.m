function [Out,max_dist] = getCellDistances(DF)
% this function produces a matrix of all cell distances between neurons 

Out = {};
max_dist  = 0;

CellIDs = DF.CellID;
active = DF.active{:,2};
Expt_list = DF.experiment_list;

for expt_idx =1:length(CellIDs)

microns_px = CellIDs{expt_idx}.micronsPerPixelX;
cell_id_microns = CellIDs{expt_idx}.ptsIdx(:,2:3) * microns_px;

active_idx = active(Expt_list == expt_idx)> 0 ;

cell_id_microns = cell_id_microns(active_idx,:);

dist = squareform(pdist(cell_id_microns));

max_dist_animal = nanmax(dist(:));

max_dist = max( [max_dist,max_dist_animal]);
 
dist(dist == 0) = nan;

Out{expt_idx} = dist;

end