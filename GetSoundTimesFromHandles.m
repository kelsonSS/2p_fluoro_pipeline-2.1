function SoundTimes = GetSoundTimesFromHandles(handles)
% This function takes a handles object and returns the sound onset and ofset
% times 

SoundTimes = struct();


Tone_on = handles.PreStimSilence * handles.pfs +1;
Tone_off = Tone_on -1 + handles.PrimaryDuration * handles.pfs +1;

if length(handles.BackgroundNoise) == 3
    Noise_on = handles.BackgroundNoise(2) * handles.pfs +1;
    Noise_off = handles.BackgroundNoise(3) * handles.pfs +1;

else 
    Noise_on = nan;
    Noise_off = nan;
end 

SoundTimes.Tone_on = Tone_on;
SoundTimes.Tone_off = Tone_off;
SoundTimes.Noise_on = Noise_on;
SoundTimes.Noise_off = Noise_off;
SoundTimes.FPS = handles.pfs;



