function [classes,active] = SortToneNoise(Passive,DFF)
% this function will take the passive object and cluster them into tone and
% noise responding neurons and perform activity analysis on them 


%init
if~exist('DFF','var')
    DFF = Passive.DFF;
end 

Class_IDS = {'Noise','Tones','Offset'}


% activity analysis 
neurons = size(DFF,3); 
fps =  Passive.handles.pfs;
soundon = Passive.handles.PreStimSilence * fps; 
soundoff = Passive.handles.BackgroundNoise(3) * fps;


active = zeros(neurons,1); 
for cell = 1:neurons
    [~,p]=ttest2(nanmean(DFF(1:soundon,:,cell)),nanmean(DFF(soundon+1:soundoff,:,cell)));
    active(cell) = p;
end 

  % find significantly active neurons  number indicates the number of
   % stars that would be added to the figure. Zero indicates inactivity
   % while anything above 1 is active
   %
   active(active <.001) = 3;  
   active(active <.01)  = 2;  
   active(active <=.05) = 1;  
   active(active < 1  ) = 0;     
   active(isnan(active)) = 0;




% cluster analysis 
DFF2 = squeeze(mean(DFF,2));

[~, onset] = max(DFF2);
 onset(isnan(onset)) = 0 ;

% relevant time-points
tone_on = Passive.handles.PreStimSilence * fps; 
tone_off = tone_on + Passive.handles.PrimaryDuration * fps;
noise_on = Passive.handles.BackgroundNoise(2) * fps;

% clustering
classes = zeros(size(DFF2,2),1); 
classes(onset>=noise_on & onset<tone_on) = 1;
classes(onset>=tone_on & onset<tone_off) = 2 ;
classes(onset>=tone_off & onset< tone_off + fps )= 3;
% remove nonresponsive neurons 
%classes(~active) = 0 ;

   
   




end 