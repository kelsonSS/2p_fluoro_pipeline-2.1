DFF_diff = diff(Passive.DFF);
fast_sig = zeros(921,4);
bl_idx = [1 25]; % baseline frames idx
tm_idx = [29,59,89,119,149];  % noise-on , tone-on,
                                    % tone-off , noise-off,
                                    % trial-end

Baseline =  squeeze(mean(DFF_diff(bl_idx(1):bl_idx(2),:,:)));

Noise_on  =  squeeze(mean(DFF_diff(tm_idx(1):tm_idx(2),:,:)));

Tone_on  =  squeeze(mean(DFF_diff(tm_idx(2):tm_idx(3),:,:)));

Tone_off  =  squeeze(mean(DFF_diff(tm_idx(3):tm_idx(4),:,:)));

Noise_off  =  squeeze(mean(DFF_diff(tm_idx(4):tm_idx(5),:,:)));

Times = {Noise_on,Tone_on,Tone_off,Noise_off};
 for neuron= 1:size(Baseline,2)
     for timePoint = 1:4 
     fast_sig(neuron,timePoint) = ttest2(Baseline(:,neuron),...
                      Times{timePoint}(:,neuron));
     end 
 end 