function TimingInfo = getTimingInfo(ThorTimingFile, PsignalTimingFile, fps)
%This function finds the frames that correspond to the end of an image
%sequence, i.e. the end of a psignal trial
TimingInfo=[];

%Load Psignal timing file
PsignalData = load(PsignalTimingFile);
%Trial duration
Primary = get(PsignalData.exptparams.TrialObject,'PrimaryHandle');
% SoundObject=Primary.SoundObject;
if strcmp(get(PsignalData.exptparams.TrialObject,'PrimaryClass'), 'Ripples')
    Prestim = 1;
    Poststim= 2;
    Duration = 3;
    TrialDur =  Prestim+Duration+Poststim;

else
    
    try
        Prestim = get(Primary,'PreStimSilence');
        Duration = get(Primary,'Duration');
        Poststim = get(Primary,'PostStimSilence');
       
        TrialDur =  Prestim + Duration + Poststim;
    catch 
         Prestim = get(Primary.SoundObject,'PreStimSilence');
        Poststim = get(Primary.SoundObject,'PostStimSilence');
        Duration =     Primary.Duration;
        TrialDur =  Prestim+Duration+Poststim;
    end
    try 
        if strcmpi(PsignalData.exptparams.BehaveObjectClass,'Passive')
            minITI = min(get(PsignalData.exptparams.BehaveObject,'ITI'));
        else
            minITI = min(get(PsignalData.exptparams.BehaveObject,'ITIs'));
        end
    catch 
        minITI = 3
    end
end

%Load ThorImage timing file--time stamps for every frame acquired,
%big jumps at ITIs. Variance in frame timing is in the writing of the frame, not acquisition.
fileID = fopen(ThorTimingFile);
formatSpec = '%f';
%Timing of frames across entire experiment
TimingVals = fscanf(fileID,formatSpec);
fclose(fileID);
%SeqEndValsTmp indexes the end of each Sequence - find peaks that
%correspond to ITIs. Because of jitter in the timing values, the peaks in
%diff(TimingVals) do not always correspond to the expected ITIs. Thus we
%have to use a method that accounts for the probablistic jitter. Here, we
%assume that the jitter does not go below 2 seconds from the minimum ITI.


jitterfactor = 0;
while minITI-jitterfactor
    SeqEndValsTmp= findpeaks(diff(TimingVals),minITI-jitterfactor);
    %Findpeaks from Chronux Toolbox--overwrites Matlab findpeaks; length(SeqEndVals) is length(trial)-1,
    %since first value indicates end of first trial
    TimingInfo.SeqEndVals = [SeqEndValsTmp.loc; length(TimingVals)];
    %Check if #trials esitmated from imaging data corresponds to actual #trials
    TotalTrials = PsignalData.exptevents(end).Trial;
    if length(TimingInfo.SeqEndVals) ~= TotalTrials
        jitterfactor = jitterfactor + .5;
    else
        break
    end
end

if length(TimingInfo.SeqEndVals) ~= TotalTrials
    
    bb = strsplit(ThorTimingFile,'/');
    warndlg(sprintf('# Trials estimated from %s is incorrect!!! \n Check SeqEndVals jitterfactor',...
        fileparts( PsignalTimingFile) ));
end
%Need to set this bc different stim sequences can have different # of frames doesn't Assume anything
TimingInfo.tarFnums = (TrialDur*fps);
%Find frame timing for each trial; length(SeqEndVals)==#trials
FrameIdx=[];
for i = 1:length(TimingInfo.SeqEndVals)
    if i == 1
        FrameIdx = [FrameIdx; [1 min([TimingInfo.tarFnums TimingInfo.SeqEndVals(1)])]];
    else
        FrameIdx = [FrameIdx; [TimingInfo.SeqEndVals(i-1)+1 TimingInfo.SeqEndVals(i-1)+ ...
            min([TimingInfo.tarFnums TimingInfo.SeqEndVals(i)])]];
    end
end
% sometimes last frame is greater than actual number of frames
if length(TimingVals)<FrameIdx(end)
    FrameIdx(end)=length(TimingVals)
else
end

TimingInfo.FrameTiming   = TimingVals(2:end); % first timingVal is always 0
TimingInfo.FrameIdx = FrameIdx;
%TimingInfo.minITI = minITI;