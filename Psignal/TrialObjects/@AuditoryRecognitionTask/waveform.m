function [TrialSound, events, o] = waveform (o,idx)
global globalparams
o=ObjUpdate(o);
TrialIndex = get(o,'TrialIndices');
%check if probe is non-target
nontargetprobe=0;
if strcmpi(get(o,'ProbeClass'), 'Sounds')
    nontargetprobe = get(get(o,'ProbeHandle'),'NonTarget');
end
BackgroundNoise=[];
if TrialIndex(idx,2)==-1 || (TrialIndex(idx,2)==0 && nontargetprobe) %Flag for probe stimulus
    StimHandle = get(o,'ProbeHandle');
    if  strcmpi(lower(deblank(get(StimHandle,'descriptor'))),'sounds')
        Type=deblank(get(StimHandle,'Type'));
    else
        Type=[];
    end
    if strcmpi(lower(deblank(get(StimHandle,'descriptor'))),'ripples')
        [TrialSound,events]=waveform(StimHandle,TrialIndex(idx,1));
    else
        [TrialSound,events]=waveform(StimHandle,TrialIndex(idx,:));
    end
    %Background noise
    if sum(get(StimHandle,'BackgroundNoise'))~=0
        BackgroundNoise=waveform(StimHandle,-1);
        NoiseParams=get(StimHandle,'BackgroundNoise');
    end
else
    StimHandle = get(o,'PrimaryHandle'); % get Primary object
    StimHandle = set(StimHandle,'Duration',TrialIndex(idx,4));
    Type=deblank(get(StimHandle,'Type'));
    
    [TrialSound,events]=waveform(StimHandle,TrialIndex(idx,:));
    
    %Background noise
    if sum(get(StimHandle,'BackgroundNoise'))~=0
        BackgroundNoise=waveform(StimHandle,-1);
        NoiseParams=get(StimHandle,'BackgroundNoise');
    end
end

%Assign behavioral meaning
ProbeClass = get(o,'ProbeClass');
PrimaryHandle = get(o,'PrimaryHandle');
switch get(o,'PrimaryClass')
    case 'Sounds'
        Discrim = get(PrimaryHandle,'NumFreqRange')>1;
    case 'SpectrumShift'
        Discrim = get(PrimaryHandle,'NumShifts')>1;
end
if TrialIndex(idx,2)==0 %NonTarget Stimulus
    tag='NonTarget';
elseif  TrialIndex(idx,2)==-1 && strcmpi(ProbeClass,'Silence') && ~Discrim
    tag='NonTarget';
elseif TrialIndex(idx,2)==-1 && strcmp(ProbeClass,'Ripples')
    tag='NonTarget';
elseif TrialIndex(idx,2)==-1 && strcmpi(ProbeClass,'Silence') && Discrim
    tag='Probe';
elseif TrialIndex(idx,2)==-1 && ~strcmp(ProbeClass,'Silence') %Probe Stimulus
    tag='Probe';
elseif TrialIndex(idx,2)>=1  %Target Stimulus
    tag='Target';
end
events=updateEvs(events,tag);
if ~isempty(strfind(Type,'PitchStim'))
    switch TrialIndex(idx,5)
        case 1
            pitchstim='IRN';
        case 2
            pitchstim='clicktrain';
        case 3
            pitchstim='AMNoise';
        case 4
            pitchstim = 'HarmStack';
        case 5
            pitchstim = 'MissFund';
    end
    events=updateEvs(events,pitchstim);
    fprintf('Trial Type: %s \n',[tag ' ' pitchstim]);
elseif ~isempty(strfind(Type,'WhiteNoise'))
    events=updateEvs(events,'WhiteNoise');
    fprintf('Trial Type: %s \n',[tag ' ' 'WhiteNoise']);
elseif sum(get(StimHandle,'BackgroundNoise'))~=0
    SNR = ['SNR=' num2str(TrialIndex(idx,3)) 'dB'];
    events=updateEvs(events,SNR);
    bb = strsep(events(2).Note,',');
    if isobject(bb{2})
        bb{2} = get(bb{2},'descriptor');
        fprintf('Trial Type: %s \n',[tag ', ' num2str(bb{2}) ', ' bb{4}]);
    else
        fprintf('Trial Type: %s \n',[tag ', ' num2str(bb{2}) ', ' bb{4}]);
    end
else
    bb = strsep(events(2).Note,',');
    if isobject(bb{2})
        bb{2} = get(bb{2},'descriptor');
        fprintf('Trial Type: %s \n',[tag ', ' num2str(bb{2}) 'Hz']);
    else
        fprintf('Trial Type: %s \n',[tag ', ' num2str(bb{2}) 'Hz']);
    end
end

%Attenuations
switch TrialIndex(idx,2)
    case {0, -1} %NonTarget amplitude scaling
        NonTargetAttendB=get(o,'NonTargetAttenDB');
        AttendB=TrialIndex(idx,3)+NonTargetAttendB;
        AttendB=10^(AttendB/20);
        TrialSound=TrialSound*AttendB;
        tag=num2str(AttendB);
    case 1 %Target amplitude scaling
        AttendB=TrialIndex(idx,3);
        AttendB=10^(AttendB/20);
        TrialSound=TrialSound*AttendB;
        tag=num2str(AttendB);
end
events=updateEvs(events,tag);
if ~isempty(BackgroundNoise)
    AttendB=NoiseParams(1);
    AttendB=10^(AttendB/20);
    BackgroundNoise=BackgroundNoise*AttendB;
    TrialSound = TrialSound+BackgroundNoise;
end

function ev=updateEvs(ev,tag)
for i=1:length(ev)
    ev(i).Note = deblank([ev(i).Note ' ,' tag]);
end
