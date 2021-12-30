function PlotLickDependentAverage(TNBehavior_active)


data = []
for expt_idx = 1:length(TNBehavior_active)
     
    expt = TNBehavior_active{expt_idx};
    if isempty(expt.LickTiming)
        continue
    end 
    data = getDataExpt(expt,data);
    

end 

end 

function data = getDataExpt(e,data)

  handles = e.handles{1};
  early_idx = logical(handles.Early);
  lick_times = e.LickTiming.TrialLickFrames
  %get early trials  
  if any(early_idx)
      
    fr_frames =cellfun(@getFirstResponse,lick_times
    Early_fluoro = e.DFF(:,early_idx,:);
    
    
  end
    
  
    
   
  
  

end 

function r = getFirstResponse(x)

if isempty(x)
    r = nan;
else 
    r = x(1);
end 
