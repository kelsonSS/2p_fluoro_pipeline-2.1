function Out = MungeBehaviorGroupData(behavior)
% quick function to create lick latencies for an entire group
% behavior is a struct with each entry containing named entries as outputs
% from BehavioralAnalysisByAnimal
%  ie Behavior.animal_01.resultsAll 
% animal name is not critical but resultsAll and its relative location is 


SNRs = [20 10 0 -10];  
AnimalIDs  =  fieldnames(behavior);



for curr_snr = 1:length(SNRs)
	SNRLatency = [] ;
    HitRateMean = [] ;
    NumTrials = [] ;
    for expt_idx = 1:length(AnimalIDs)
        Animal_ID = AnimalIDs{expt_idx};
        if strcmp(Animal_ID,'group')
            continue
        end 
    expt = behavior.(Animal_ID).resultsAll.Combined;
    
    snr_idx = find(expt.SNR == SNRs(curr_snr));
    if isempty(snr_idx)
        continue
    end 
     
    
    SNRLatency= cat(1,SNRLatency, expt.LickLatency{snr_idx});
    HitRateMean = cat(1,HitRateMean,expt.HitRateMean(snr_idx));
    NumTrials = cat(1,NumTrials,expt.NumTrials(snr_idx));
    
    end 
    % packing
    snr_name = sprintf('SNR_%d',SNRs(curr_snr));
    snr_name = strrep(snr_name,'-','minus_');
    Out.(snr_name).SNRLatency = SNRLatency;
    Out.(snr_name).HitRateMean = HitRateMean;
    Out.(snr_name).NumTrials = NumTrials;
    
end 
    


