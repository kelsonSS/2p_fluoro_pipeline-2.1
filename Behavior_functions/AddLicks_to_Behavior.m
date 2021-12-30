function TNBehavior = AddLicks_to_Behavior(TNBehavior,Licks)


for expt_idx = 1:length(TNBehavior)
        
        L = Licks{expt_idx};
        
        if isempty(L)
          TNBehavior{expt_idx}.LickTiming = {};
        else 
            
        try
        TNBehavior{expt_idx}.LickTiming = L.Out; 
        catch
         TNBehavior{expt_idx}.LickTiming = L;
        end 
        end 
end

