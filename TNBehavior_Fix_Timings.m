function TNBehavior = TNBehavior_Fix_Timings(TNBehavior)
% this function takes the prestim silence and trims it to one second 

for expt_idx = 1:length(TNBehavior)
    if TNBehavior{expt_idx}.handles{1}.PreStimSilence >1
        TNBehavior{expt_idx}.DFF     = TNBehavior{expt_idx}.DFF(31:end,:,:);
        TNBehavior{expt_idx}.DFF_Z   = TNBehavior{expt_idx}.DFF_Z(31:end,:,:);
        TNBehavior{expt_idx}.DFF_norm = TNBehavior{expt_idx}.DFF_norm(31:end,:,:);
        TNBehavior{expt_idx}.handles{1}.PreStimSilence  =...
            TNBehavior{expt_idx}.handles{1}.PreStimSilence -1; 
        
    end 
end