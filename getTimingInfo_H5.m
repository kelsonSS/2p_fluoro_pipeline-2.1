function TimingInfo = getTimingInfo_H5(ThorSyncFile, PsignalTimingFile, fps,xml)
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
try
gate=h5read(ThorSyncFile,'/AI/ai4'); % trial gate signal
fo=h5read(ThorSyncFile,'/DI/Frame Out');
catch    
    return
end 

try
fc=h5read(ThorSyncFile,'/CI/Frame Counter'); % may be unreliable
catch
    fc = fo;
end 
    
nchannels = xml.format{2}(3);

%first_frame = find(fc,1);
    

% find onsets and offsets by looking at the derivative of the gate 
on  = findpeaks(diff(gate),1);
off = findpeaks(diff(gate * -1),1);

on  = floor(on.loc/1000);
off = floor(off.loc/1000);




% readjusts frame timing if there was a delay before the first frame 


firstframe = floor(find(fc,1)/1000);
on = on - firstframe + 1;
off = off - firstframe + 1;

if on(1) == 0
    on(1) = 1; 
    off = off+1;
end 


%sanity checks#
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

 TotalTrials = PsignalData.exptevents(end).Trial ; %

%Find frame timing for each trial; length(SeqEndVals)==#trials
FrameIdx(:,1) = on;
FrameIdx(:,2) = on + TimingInfo.tarFnums;

TimingInfo.FrameIdx = FrameIdx;
TimingInfo.SeqEndVals = FrameIdx(:,2);
TimingInfo.FrameTiming = frameseconds;
%TimingInfo.ITI = ITI;


