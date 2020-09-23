function out = CreateBehaviorCellList(Behavior)

if size(Behavior,2) >1
    out.Passive = CreateBehaviorCellList(Behavior(:,1));
    out.Active = CreateBehaviorCellList(Behavior(:,2));
    return 
end 

out = struct();
out.Clean_idx = [];
out.active = [];
out.Cells_used = {};
out.NCells = [];
for ii = 1:length(Behavior)
    expt = Behavior{ii};
    active_idx = expt.active{:,2} > 0;
    out.Cells_used{ii} = find(expt.Clean_idx & active_idx);
    out.NCells = cat(1,out.NCells, length(out.Cells_used{ii}) );
    % append fields to master list
    out.Clean_idx = cat(1,out.Clean_idx,expt.Clean_idx);
    out.active = cat(1,out.active,expt.active);
  
end

        

