function close_idx = CompareBestFrequencyToBehaviorTarget(TNbehavior)



for expt = 1:size(TNbehavior,1);
    
   responsive =  TNbehavior{expt,1}.active.Activity > 0 ; 
   BFs = TNbehavior{expt,1}.BestFrequencies{1}(responsive);
   
   active_handles = TNbehavior{expt,2}.handles{1};
   
   active_freq = active_handles.uFreqs;
    
   % find all Freqs within 1 octave of BF
    close_idx{expt,1} = abs(log2(BFs)-log2(active_freq)) < 1;
    
end