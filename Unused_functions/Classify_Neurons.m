function out = Classify_Neurons(DFF,handles)

% this function takes in a matrix in the form 
%      
%        Frame X Trial X Neuron 
% and classifies each neuron into how it responded to given stimuli 

% Average trials if individual trial data given
  
 
 % Settings currently hard coded here 
 silence_end = 30;  
 noise_on =  31;
 tone_on = 61 ;
 tone_off = 90;
 noise_off = 120; 
                                     
 % next we will create the average response of  each neuron
  baseline_mu = squeeze(nanmean(DFF(1:noise_on-1,: , :) )) ;
  noise_mu = squeeze(nanmean( DFF(noise_on:tone_on,: , :) )) ;
  tone_mu = squeeze(nanmean(DFF(tone_on:tone_off,:,:) )) ;
  off_mu  = squeeze(nanmean(DFF(tone_off:noise_off,:,:) ));
 
 for neuron = 1:size(DFF,3)
      ii = neuron;
     
     Noise(ii) = ttestsig(baseline_mu(:,neuron),noise_mu(:,neuron)); 
     Tone(ii)  = ttestsig(baseline_mu(:,neuron),tone_mu(:,neuron));
     TonePreferring(ii) = ttestsig(noise_mu(:,neuron),tone_mu(:,neuron));
    % ToneAlone(ii)  = 
    % NoiseAlone(ii) =
     NegativeNoise(ii) = min(nanmean(DFF(noise_on:tone_on,:,neuron),2)) < 0;  
     NegativeTone(ii) = min(nanmean(DFF(tone_on:tone_off,:,neuron),2)) < 0;  
     
 end 
 
 Sig.Noise  = Noise;
 Sig.Tone = Tone;
 Sig.TonePreferring = TonePreferring;
 Sig.NegativeNoise = NegativeNoise;
 Sig.NegativeTone = NegativeTone;
 
 
end 

function bool = ttestsig(s1,s2,sig)
if ~exist('sig', 'var')
    sig =.01;
end 
 [~,p ] = ttest2(s1,s2); 
  
 bool = [ p <= sig ];
 
end 









