%young_oldPassiveFig2

% get indicies of sounds played at indicated SNR 
tones_idx = FRA_young_noise.FreqLevelOrder{:,2} == inf;
noise20SNR_idx =FRA_young_noise.FreqLevelOrder{:,2} ==70 ;

CDFStats = {};
stats_idx = 1;
% 
Types = {'mean','max'};
Levels = {'Tones','+20db SNR'};
level_idx ={tones_idx,noise20SNR_idx};

for type = 1:2
    for lvl = 1:2
        % get CDF 
      young =  PlotFluoroCDF(FRA_young_noise,Types{type},level_idx{lvl});
        old =  PlotFluoroCDF(FRA_old_noise,Types{type},level_idx{lvl});
        % get stats using KS-test2 
        CDFStats{stats_idx,1} = sprintf( '%s , %s', Types{type},Levels{lvl});
        [~,CDFStats{stats_idx,2}] = kstest2(young,old);
        stats_idx = stats_idx+1;
    end 
end 

clear Types 
clear Levels
clear level_idx  
clear stats_idx
clear tones_idx
clear noise20SNR_idx

