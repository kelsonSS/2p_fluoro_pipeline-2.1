function GeneratePsignalMatrix(input)
%This script generates Psignal matrices, P, which are Time-index x Trial x M,
%where M is a task-related feature. P(i,t,m) is equal to eqither 1, 0, time value
%in seconds, or feature value, with 1 indicating  membership of a given feature
%or time value of event. Lick have a value of 1 whenever a lick occured.
fps = input.expectedFPS;
for e = 1:length(input.expname)
    PsignalFile = input.psignalfiles{e};
    bb = strsplit(input.path,'\');
    InPath   = fullfile(input.path , input.expname{e}, input.psignalfiles{e}) 
    SavePath = input.savepath;
    SavePath = [SavePath '\'  bb{5}  '\' input.expname{e} '\'];
    copyfile ( InPath, SavePath);
    PsignalFile=[SavePath PsignalFile];
    PsignalData = load(PsignalFile);
    Primary = get(PsignalData.exptparams.TrialObject,'PrimaryHandle');
    try 
        SoundObject = get(Primary, 'SoundObject')
    catch    
    SoundObject = Primary.SoundObject
    end
    globalparams = PsignalData.globalparams;
    disp(['Processing ' PsignalFile])
    TotalTrials = PsignalData.exptevents(end).Trial;
    TrialLen = PsignalData.exptevents(end).StopTime;
    TrialObject =PsignalData.exptparams.TrialObject;
   try
    OveralldB = TrialObject.OveralldB;
   catch
    OveralldB = get(TrialObject, 'OveralldB');
   end 
    try 
        StimOnset = ceil(get(SoundObject,'PreStimSilence')*fps)+1;
    catch
        StimOnset = ceil(getfield(SoundObject,'PreStimSilence')*fps)+1;
    end
    try
        StimOffset = (get(SoundObject,'PreStimSilence')+get(Primary,'Duration'))*fps;
    catch
        StimOffset = (getfield(SoundObject,'PreStimSilence')+getfield(Primary,'Duration'))*fps;
    end
    
    %Trialobject dependent features
    if ~isempty(strfind(PsignalFile,'SHP'))
        BehaviorObject = get(PsignalData.exptparams.BehaveObject);
        %Features
        PsignalMatrix=[];
        for i = 1:7
            PsignalMatrix.TagNames{i,1} = i;
        end
        PsignalMatrix.TagNames{1,2} = 'Hit';
        PsignalMatrix.TagNames{2,2} = 'Miss';
        PsignalMatrix.TagNames{3,2} = 'FirstResponse';
        PsignalMatrix.TagNames{4,2} = 'TargetLick';
        PsignalMatrix.TagNames{5,2} = 'StimOnset';
        PsignalMatrix.TagNames{6,2} = 'StimOffset';
        PsignalMatrix.TagNames{7,2} = 'StimLevel';
        %Preallocation of P
        PsignalMatrix.Tags = zeros(floor(fps*TrialLen),TotalTrials,size(PsignalMatrix.TagNames(:,2),1));
        %Performance
        if isfield(PsignalData.exptparams,'Performance')
            for i = 1:PsignalData.exptevents(end).Trial
                switch PsignalData.exptparams.Performance(i).ThisTrial
                    case 'Hit'
                        idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'Hit'));
                        PsignalMatrix.Tags(StimOnset,i,idx) = 1;
                    case 'Miss'
                        idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'Miss'));
                        PsignalMatrix.Tags(StimOnset,i,idx) = 1;
                end
            end
            %First response
            for i = 1:PsignalData.exptevents(end).Trial
                idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'FirstResponse'));
                if ~isnan(PsignalData.exptparams.FirstResponse(i,1))
                    PsignalMatrix.Tags(max([1 ceil(PsignalData.exptparams.FirstResponse(i,1).*fps)]),i,idx) = 1;
                end
            end
            %Licking
            [t,trial,Note,toff,StimIndex] = evtimes(PsignalData.exptevents,'Stim*');
            Ti=0;
            Pi=0;
            NTi=0;
            for i = 1:PsignalData.exptevents(end).Trial
                bb = strsplit(PsignalData.exptevents(StimIndex(i)).Note,',');
                StimType = bb{3}(1:end);
                switch StimType
                    case 'Target'
                        Ti = Ti+1;
                        idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'Target'));
                        PsignalMatrix.Tags(StimOnset,i,idx) = 1;
                        ResponseData = PsignalData.exptparams.AverageResponse.tar(Ti,:);
                        Lfs = globalparams.HWparams.fsAI;
                        t=0:1/Lfs:(length(ResponseData)/Lfs)-(1/Lfs);
                        ResponseFrames=max(1,round(t(find(ResponseData)).*fps));
                        idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'TargetLick'));
                        PsignalMatrix.Tags(ResponseFrames,i,idx) = 1;
                end
            end
        end
        %Stimulus
        [t,trial,Note,toff,StimIndex] = evtimes(PsignalData.exptevents,'Stim*');
        for i = 1:PsignalData.exptevents(end).Trial
            bb = strsplit(PsignalData.exptevents(StimIndex(i)).Note,',');
            if isempty(strfind(Note{i},'Silence'))
                idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'StimOnset'));
                PsignalMatrix.Tags(StimOnset,i,idx) = 1;
                idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'StimOffset'));
                PsignalMatrix.Tags(StimOffset,i,idx) = 1;
            end
        end
     elseif ~isempty(strfind(PsignalFile,'ART'))
         BehaviorObject = get(PsignalData.exptparams.BehaveObject);
         BehaviorClass  = PsignalData.exptparams.TrialObjectClass;
         if isempty(strmatch(BehaviorClass, 'Passive'))
         EarlyWindowOnset = StimOnset+floor(BehaviorObject.EarlyWindow(1)*fps);
         EarlyWindowOffset = StimOnset+floor(BehaviorObject.EarlyWindow(2)*fps)-1;
         ResponseWindowOnset = EarlyWindowOffset+1;
         ResponseWindowOffset =  ResponseWindowOnset+diff(BehaviorObject.ResponseWindow)*fps-1;
         %Features
         PsignalMatrix=[];
         for i = 1:20
             PsignalMatrix.TagNames{i,1} = i;
         end
         PsignalMatrix.TagNames{1,2} = 'Hit';
         PsignalMatrix.TagNames{2,2} = 'Miss';
         PsignalMatrix.TagNames{3,2} = 'FalseAlarm';
         PsignalMatrix.TagNames{4,2} = 'CorrectRejection';
         PsignalMatrix.TagNames{5,2} = 'Early';
         PsignalMatrix.TagNames{6,2} = 'FirstResponse';
         PsignalMatrix.TagNames{7,2} = 'TargetLick';
         PsignalMatrix.TagNames{8,2} = 'NonTargetLick';
         PsignalMatrix.TagNames{9,2} = 'ProbeLick';
         PsignalMatrix.TagNames{10,2} = 'StimOnset';
         PsignalMatrix.TagNames{11,2} = 'StimOffset';
         PsignalMatrix.TagNames{12,2} = 'StimFrequency';
         PsignalMatrix.TagNames{13,2} = 'StimLevel';
         PsignalMatrix.TagNames{14,2} = 'Target';
         PsignalMatrix.TagNames{15,2} = 'NonTarget';
         PsignalMatrix.TagNames{16,2} = 'Probe';
         PsignalMatrix.TagNames{17,2} = 'EarlyWindowOnset';
         PsignalMatrix.TagNames{18,2} = 'EarlyWindowOffset';
         PsignalMatrix.TagNames{19,2} = 'ResponseWindowOnset';
         PsignalMatrix.TagNames{20,2} = 'ResponseWindowOffset';
         %Preallocation of P
         PsignalMatrix.Tags = zeros(floor(fps*TrialLen),TotalTrials,size(PsignalMatrix.TagNames(:,2),1));
         %Performance
         if isfield(PsignalData.exptparams,'Performance')
             for i = 1:PsignalData.exptevents(end).Trial
                 switch PsignalData.exptparams.Performance(i).ThisTrial
                     case 'Hit'
                         idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'Hit'));
                         PsignalMatrix.Tags(StimOnset,i,idx) = 1;
                     case 'Miss'
                         idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'Miss'));
                         PsignalMatrix.Tags(StimOnset,i,idx) = 1;
                     case 'FalseAlarm'
                         idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'FalseAlarm'));
                         PsignalMatrix.Tags(StimOnset,i,idx) = 1;
                     case 'CorrectReject'
                         idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'CorrectReject'));
                         PsignalMatrix.Tags(StimOnset,i,idx) = 1;
                     case 'Early'
                         idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'Early'));
                         PsignalMatrix.Tags(StimOnset,i,idx) = 1;
                 end
             end
             %First response
             for i = 1:PsignalData.exptevents(end).Trial
                 idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'FirstResponse'));
                 if ~isnan(PsignalData.exptparams.FirstResponse(i,1))
                     PsignalMatrix.Tags(max([1 ceil(PsignalData.exptparams.FirstResponse(i,1).*fps)]),i,idx) = 1;
                 end
             end
             %Licking
             [t,trial,Note,toff,StimIndex] = evtimes(PsignalData.exptevents,'Stim*');
             Ti=0;
             Pi=0;
             NTi=0;
             for i = 1:PsignalData.exptevents(end).Trial
                 bb = strsplit(PsignalData.exptevents(StimIndex(i)).Note,',');
                 StimType = bb{3}(1:end-1);
                 switch StimType
                     case 'Target'
                         Ti = Ti+1;
                         idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'Target'));
                         PsignalMatrix.Tags(StimOnset,i,idx) = 1;
                         ResponseData = PsignalData.exptparams.AverageResponse.tar(Ti,:);
                         Lfs = globalparams.HWparams.fsAI;
                         t=0:1/Lfs:(length(ResponseData)/Lfs)-(1/Lfs);
                         ResponseFrames=max(1,round(t(find(ResponseData)).*fps));
                         idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'TargetLick'));
                         PsignalMatrix.Tags(ResponseFrames,i,idx) = 1;
                     case 'NonTarget'
                         NTi = NTi+1;
                         idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'NonTarget'));
                         PsignalMatrix.Tags(StimOnset,i,idx) = 1;
                         ResponseData = PsignalData.exptparams.AverageResponse.nontar(NTi,:);
                         Lfs = globalparams.HWparams.fsAI;
                         t=0:1/Lfs:(length(ResponseData)/Lfs)-(1/Lfs);
                         ResponseFrames=max(1,round(t(find(ResponseData)).*fps));
                         idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'NonTargetLick'));
                         PsignalMatrix.Tags(ResponseFrames,i,idx) = 1;
                     case 'Probe'
                         Pi = Pi+1;
                         idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'Probe'));
                         PsignalMatrix.Tags(StimOnset,i,idx) = 1;
                         ResponseData = PsignalData.exptparams.AverageResponse.probe(Pi,:);
                         Lfs = globalparams.HWparams.fsAI;
                         t=0:1/Lfs:(length(ResponseData)/Lfs)-(1/Lfs);
                         ResponseFrames=max(1,round(t(find(ResponseData)).*fps));
                         idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'ProbeLick'));
                         PsignalMatrix.Tags(ResponseFrames,i,idx) = 1;
                 end
             end
             %Behavior Windows
             for i = 1:PsignalData.exptevents(end).Trial
                 idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'EarlyWindowOnset'));
                 PsignalMatrix.Tags(EarlyWindowOnset,i,idx) = 1;
                 idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'EarlyWindowOffset'));
                 PsignalMatrix.Tags(EarlyWindowOffset,i,idx) = 1;
                 idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'ResponseWindowOnset'));
                 PsignalMatrix.Tags(ResponseWindowOnset,i,idx) = 1;
                 idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'ResponseWindowOffset'));
                 PsignalMatrix.Tags(ResponseWindowOffset,i,idx) = 1;
             end
         end
         %Stimulus
         [t,trial,Note,toff,StimIndex] = evtimes(PsignalData.exptevents,'Stim*');
         for i = 1:PsignalData.exptevents(end).Trial
             bb = strsplit(PsignalData.exptevents(StimIndex(i)).Note,',');
             if isempty(strfind(Note{i},'Silence'))
                 idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'StimOnset'));
                 PsignalMatrix.Tags(StimOnset,i,idx) = 1;
                 idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'StimOffset'));
                 PsignalMatrix.Tags(StimOffset,i,idx) = 1;
                 idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'StimFrequency'));
                 PsignalMatrix.Tags(StimOnset,i,idx) = bb{2};
                 idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'StimLevel'));
                 PsignalMatrix.Tags(StimOnset,i,idx) = OveralldB;
             end
         end

       %Features
       PsignalMatrix=[];
       for i = 1:4
           PsignalMatrix.TagNames{i,1} = i;
       end
       PsignalMatrix.TagNames{1,2} = 'StimOnset';
       PsignalMatrix.TagNames{2,2} = 'StimOffset';
       PsignalMatrix.TagNames{3,2} = 'StimFrequency';
       PsignalMatrix.TagNames{4,2} = 'StimLevel';
       %Preallocation of P
       PsignalMatrix.Tags = zeros(floor(fps*TrialLen),TotalTrials,size(PsignalMatrix.TagNames(:,2),1));
       %Stimulus
       [t,trial,Note,toff,StimIndex] = evtimes(PsignalData.exptevents,'Stim*');
       for i = 1:TotalTrials
           bb = strsplit(PsignalData.exptevents(StimIndex(i)).Note,',');
           idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'StimOnset'));
           PsignalMatrix.Tags(StimOnset,i,idx) = 1;
           idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'StimOffset'));
           PsignalMatrix.Tags(StimOffset,i,idx) = 1;
           idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'StimFrequency'));
           PsignalMatrix.Tags(StimOnset,i,idx) = str2num(bb{2});
           idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'StimLevel'));
           if strcmp(bb{end},'OptTrigOn')
               PsignalMatrix.Tags(StimOnset,i,idx) = OveralldB-str2num(bb{end-1}); %changed because of OptTrig option...
           else
               PsignalMatrix.Tags(StimOnset,i,idx) = OveralldB-str2num(bb{end});
           end
       end
       
       
       
       elseif  ~isempty(strfind(PsignalFile,'TORC'))
       %Features
       PsignalMatrix=[];
       for i = 1:4
           PsignalMatrix.TagNames{i,1} = i;
       end
        PsignalMatrix.TagNames{1,2} = 'StimOnset';
        PsignalMatrix.TagNames{2,2} = 'StimOffset';
        PsignalMatrix.TagNames{3,2} = 'StimFrequency';
        PsignalMatrix.TagNames{4,2} = 'StimLevel';
        %Preallocation of P
        PsignalMatrix.Tags = zeros(floor(fps*TrialLen),TotalTrials,size(PsignalMatrix.TagNames(:,2),1));
        %Stimulus
        [t,trial,Note,toff,StimIndex] = evtimes(PsignalData.exptevents,'Stim*');
        for i = 1:TotalTrials
            bb = strsplit(PsignalData.exptevents(StimIndex(i)).Note,',');
            bb = strsplit(bb{2},'_')
            idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'StimOnset'));
            PsignalMatrix.Tags(StimOnset,i,idx) = 1;
            idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'StimOffset'));
            PsignalMatrix.Tags(StimOffset,i,idx) = 1;
            idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'StimFrequency'));
            PsignalMatrix.Tags(StimOnset,i,idx) = str2num(bb{end});
            idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'StimLevel'));
            PsignalMatrix.Tags(StimOnset,i,idx) = OveralldB%-str2num(bb{end});
        end
        elseif  ~isempty(strfind(PsignalFile,'Timbre'))
        %Features
        PsignalMatrix=[];
        for i = 1:4
            PsignalMatrix.TagNames{i,1} = i;
        end
        PsignalMatrix.TagNames{1,2} = 'StimOnset';
        PsignalMatrix.TagNames{2,2} = 'StimOffset';
        PsignalMatrix.TagNames{3,2} = 'StimFrequency';
        PsignalMatrix.TagNames{4,2} = 'StimLevel';
        %Preallocation of P
        PsignalMatrix.Tags = zeros(floor(fps*TrialLen),TotalTrials,size(PsignalMatrix.TagNames(:,2),1));
        %Stimulus
        [t,trial,Note,toff,StimIndex] = evtimes(PsignalData.exptevents,'Stim*');
        for i = 1:TotalTrials
            bb = strsplit(PsignalData.exptevents(StimIndex(i)).Note,',');
            bb = strsplit(bb{2},'_')
            idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'StimOnset'));
            PsignalMatrix.Tags(StimOnset,i,idx) = 1;
            idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'StimOffset'));
            PsignalMatrix.Tags(StimOffset,i,idx) = 1;
            idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'StimFrequency'));
            PsignalMatrix.Tags(StimOnset,i,idx) = str2num(bb{end-1});
            idx = find(strcmpi(PsignalMatrix.TagNames(:,2),'StimLevel'));
            PsignalMatrix.Tags(StimOnset,i,idx) = OveralldB%-str2num(bb{end});
        end
        
    end
        save([SavePath 'PsignalMatrix.mat'],'PsignalMatrix')
    end
end