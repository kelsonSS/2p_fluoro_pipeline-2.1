function TNBehavior = TNBehavior_Fix_Timings(TNBehavior)
% this function takes the prestim silence and trims it to one second 

for expt_idx = 1:length(TNBehavior)
    
    threshold = 64000;
    above_thresh_idx = TNBehavior{expt_idx}.handles{1}.Freqs >= threshold; 
    
    if any(above_thresh_idx)
        keep_idx = ~above_thresh_idx;
        
        TNBehavior{expt_idx}.DFF     = TNBehavior{expt_idx}.DFF(:,keep_idx,:);
        TNBehavior{expt_idx}.DFF_Z   = TNBehavior{expt_idx}.DFF_Z(:,keep_idx,:);
        TNBehavior{expt_idx}.DFF_norm = TNBehavior{expt_idx}.DFF_norm(:,keep_idx,:);
        TNBehavior{expt_idx}.handles{1}.FreqLevelOrder  =...
            TNBehavior{expt_idx}.handles{1}.FreqLevelOrder(keep_idx); 
         TNBehavior{expt_idx}.handles{1}.Freqs  =...
            TNBehavior{expt_idx}.handles{1}.Freqs(keep_idx); 
         TNBehavior{expt_idx}.handles{1}.Levels  =...
            TNBehavior{expt_idx}.handles{1}.Levels(keep_idx); 
        
    end 
end