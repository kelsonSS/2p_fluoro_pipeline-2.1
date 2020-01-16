function  Corrs = Corrs_by_temporal_timing_behavior(TN)
% This script extends the functionality of Correlations function to extract
% the values for different subsets of neurons ( currently written to be
% four sets named below)

Names = {'All','Tone_on','Tone_Off','ramping'};
DFF_mu = squeeze(nanmean(TN.DFF,2));
[~,onsets]  = max(DFF_mu);

on_idx = onsets > 30 & onsets <= 60;
off_idx  = onsets > 60 & onsets <= 90;
ramping_idx = onsets > 90 & onsets<= 120;

Indicies = {logical(onsets);on_idx;off_idx;ramping_idx};
Corrs = cell(4,2);

for subgroup = 1:4    
    TN_temp = TN;
    TN_temp.DFF = TN.DFF(:,:,Indicies{subgroup});
    TN_temp.experiment_list = TN.experiment_list(Indicies{subgroup});
    TN_temp.active = TN.active(Indicies{subgroup},:);
    
    
    Corrs{subgroup,1} = Correlations(TN_temp);
    Corrs(subgroup,2) = Names(subgroup);
    
 end 
 
 clear TN_temp