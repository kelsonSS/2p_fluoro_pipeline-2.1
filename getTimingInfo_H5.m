function TimingInfo = getTimingInfo_H5(ThorSyncFile, PsignalTimingFile, fps,xml,num_frames_actual)
%This function finds the frames that correspond to trials

if ~exist('fps','var')
    fps = 30 ; 
 
end
     

%% Load Psignal timing file and Relevant Info
PsignalData = load(PsignalTimingFile);

try
    Primary = get(PsignalData.exptparams.TrialObject,'PrimaryHandle');
    SoundObject=get(Primary,'SoundObject');
catch
    Primary = PsignalData.exptparams.TrialObject.PrimaryHandle;
    SoundObject=Primary.SoundObject;
end

%% Trial Parameters
try
    Prestim = get(SoundObject,'PreStimSilence');
    Duration = get(Primary,'Duration');
    Poststim = get(SoundObject,'PostStimSilence');
    TrialDur =  Prestim+Duration+Poststim;
    TimingInfo.tarFnums = (TrialDur*fps);

    %ITI = get(PsignalData.exptparams.BehaveObject,'ITI');
catch
    Prestim= SoundObject.PreStimSilence;
    Duration = Primary.Duration;
    Poststim = SoundObject.PostStimSilence;
    TrialDur = Prestim+Duration+Poststim;
    %Need to set this bc different stim sequences can have different # of frames doesn't Assume anything
   
    %ITI = PsignalData.exptparams.BehaveObject.ITI;
end

%% Extract Timing Info from h5 file
fprintf('loading %s \n',ThorSyncFile) 
try
fo=h5read(ThorSyncFile,'/DI/Frame Out');
catch    
  fprintf('%s is bad /n', ThorSyncFile)
    return
end 

try
fc=h5read(ThorSyncFile,'/CI/Frame Counter'); % may be unreliable
catch
    fc = fo;
end 
    
try
    gate=h5read(ThorSyncFile,'/AI/ai5'); % trial gate signal
catch
    try
        gate=h5read(ThorSyncFile,'/AI/PsignalGate');
    catch
        try
        gate=h5read(ThorSyncFile,'/AI/ai4');
        catch
             fprintf('%s is bad /n', ThorSyncFile)
        end
        
    end
end



nchannels = xml.format{2}(3);

first_frame = find(fc,1);

% find onsets and offsets by looking at the derivative of the gate 

TotalTrials = PsignalData.exptevents(end).Trial ; %
peak_val = .521;

if  any(gate(50:100) > 1)  % if gate starts high fix gate trigger to show it going high  
    gate(1:10) = 0;
else 
    gate(1:100) = 0; 
end 

on =  findpeaks(diff(gate),.5);
off = findpeaks(diff(gate * -1),.5);
on = on.loc;
off = off.loc;



if length(on) ~=TotalTrials
    delta_peak = .1; % how much to change peak val by 
    pv_lower_flag= 0;  % if the last trial was too small
    pv_higher_flag = 0; % if last trial was too large
end 
while length(on) ~= TotalTrials
    m = length(on);
    g = length(gate);
    if m > 5000 || peak_val < 0.001 || peak_val > 2 || delta_peak < 5e-2
        warning('trial error manual inspection needed \n %s \n',ThorSyncFile)
        pause
        return
    end
        
    
    
    % if we have switched directions decrease the delta size
    if pv_lower_flag && pv_higher_flag
        delta_peak = delta_peak/2;
        pv_lower_flag = 0;
        pv_higher_flag = 0;
    end 
    if m < TotalTrials
        pv_lower_flag = 1;
        peak_val= peak_val - delta_peak;
        fprintf('too few estimated trials \n')   
        fprintf('estimated trials:%d actual: %d peak val: %.2f Delta: %.2f \n',...
            m,TotalTrials,peak_val,delta_peak)
    end 
    
      if m > TotalTrials
        pv_higher_flag = 1;
        peak_val= peak_val + delta_peak;
        fprintf('too many estimated trials \n')   
        fprintf('estimated trials:%d actual: %d peak val: %.2f Delta: %.2f \n',...
            m,TotalTrials,peak_val,delta_peak)
      end
    
      
    on =  findpeaks(diff(gate),peak_val);
    % ensure on peaks are acutally on

    on = on.loc;
    on_idx = gate(min(on+9,g)) - gate(max(on-10,1)) > 0;
    on = on(on_idx);
    
    off = findpeaks(diff(gate * -1),peak_val);
    %ensure off peaks are off  
    off = off.loc;
    off_idx = gate(max(off-10,1)) - gate(min(off+9,g)) > 0;
    off = off(off_idx);
    
    
end  
    



 % fix timing params to seconds    
on  = floor(on/1000);
off = floor(off/1000);    
    
if find(fc,1) == 1 
    on = cat(1,1, on);
    off = off(2:end);
end 









% readjusts frame timing if there was a delay before the first frame 


firstframe = floor(find(fc,1)/1000);
on = on - firstframe + 1;
off = off - firstframe + 1;



%% sanity checks

if num_frames_acutal < TotalTrials * TimingInfo.TarFnums -1
    warning('trial has too few frames, skipping')
    return
end 


num_frames_predicted = findpeaks(diff(fc),1);
num_frames_predicted = length(num_frames_predicted.loc);

% if more frames in actual than predicted, shift times
% to reflect true start time 
try
if num_frames_predicted < num_frames_actual
    
   shift =  num_frames_actual - num_frames_predicted;
   on =  on+shift;
   off = off+shift;
    sprintf('%s, %s' ,num_frames_predicted,num_frames_actual);
end
catch
    warning('num_frames_error')
end 
% if frames are off by one shift first start frame by 1 
% indicates one frame dropped at beginning 
if any(size(on) == 0)
    on(1) = 1;
end 

if num_frames_predicted == num_frames_actual + 1 || on(1) == 0 
        on(1) = 1;  
       % off(1) = off(1)+1;
end 
    
%  change on and off to reflect frame times if in stimulus mode 
if  strcmp(xml.Streaming,'0')
  [on,off] =  extract_h5_stimulus_mode(on,off);
end



% frames_actual = findpeaks(diff(fo(on(1)*1000:off(1)*1000 )),1);
% bad_frame_flag = length(frames_actual.loc) ~= TimingInfo.tarFnums ;
% 
% bad_timing_flag = (floor(mean(off-on)) ~= TimingInfo.tarFnums);
% 
% % sanity check
% if bad_timing_flag || bad_frame_flag
%     warning(sprintf('%s may be wrong file! manual inspection needed',ThorSyncFile))
% end 


% timing extraction

frameseconds = findpeaks(diff(fc ),1);
frameseconds = frameseconds.loc /1000/30;

%% Create TimingInfo Structure

 

 
%Find frame timing for each trial; length(SeqEndVals)==#trials
FrameIdx(:,1) = on;
FrameIdx(:,2) = on + TimingInfo.tarFnums - 1 ;

TimingInfo.FrameIdx = FrameIdx;
TimingInfo.SeqEndVals = FrameIdx(:,2);
TimingInfo.FrameTiming = frameseconds;
%TimingInfo.ITI = ITI;


function [on off] = extract_h5_stimulus_mode(on,off)



if num_frames_actual == TimingInfo.tarFnums  * TotalTrials
    on = 1 + TimingInfo.tarFnums .*(0:TotalTrials-1);
    off = on + TimingInfo.tarFnums - 1; 
return
end




for ii = 1:TotalTrials
   
    

    
    

%% Create TimingInfo Structure



 

end 

    

end

end 





        


