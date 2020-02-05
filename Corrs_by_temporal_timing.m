function  Corrs = Corrs_by_temporal_timing(TN)
% This script extends the functionality of Correlations function to extract
% the values for different subsets of neurons ( currently written to be
% four sets named below)

Names = {'All','Noise','Tone_on','Tone_off'};
DFF_mu = squeeze(nanmean(TN.DFF,2));
[~,onsets]  = max(DFF_mu);

Noise_idx = onsets > 30 & onsets <= 60;
Tone_idx  = onsets > 60 & onsets <= 90;
Offset_idx = onsets > 90 & onsets<= 120;

Indicies = {logical(onsets);Noise_idx;Tone_idx;Offset_idx};
Corrs = cell(4,2);

m = length(Indicies);
for subgroup1 = 1:4
    for subgroup2 = 1:4
        if subgroup1>subgroup2
            % upper triangular form
            continue
        end 
        group_name = sprintf('%s, %s ',Names{subgroup1},Names{subgroup2});
    figure 
    title(group_name)
    
    TN_temp = TN;
    combined_index = Indicies{subgroup1} | Indicies{subgroup2};
    TN_temp.DFF = TN.DFF(:,:,combined_index);
    TN_temp.experiment_list = TN.Experiment_list(combined_index);
    TN_temp.active = TN.Active(combined_index,:);
    
    
    Corrs{subgroup1,subgroup2} = Correlations(TN_temp);
      Corrs{subgroup1,m+1} = group_name;
    
    end
    
end 
 
 clear TN_temp