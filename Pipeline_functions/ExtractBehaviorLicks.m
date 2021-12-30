function ExtractBehaviorLicks(input,debug)

tic;



if ~exist('debug','var')
    debug= 0;
end 

    TimingInfo = struct('Errors',[]);
    %% path prep
 
    PsignalFileName = input.psignalfiles;
    SavePath = input.path;
    LocalPath = input.path;
        
    TimingFile = fullfile(LocalPath,'TimingInfo.mat');
    xmlfile = fullfile(LocalPath,'Experiment.xml');
    PsignalFile= fullfile(LocalPath,PsignalFileName);
      fps =input.expectedFPS ;
     if exist(char(fullfile(LocalPath,'redchannel.mat')),'file')
         fps = fps/2;
     end 
     
     try
     xml = get_options_from_xml(xmlfile);
     catch
       warning('bad xml file')
         return
     end 
     
    %% File Prep
    if iscell(PsignalFile)
        PsignalFile = PsignalFile{1};
    end 
    
     if ~exist(PsignalFile,'file')
        TimingInfo.Errors = LocalPath;
         return
    end
    
    
   % if file was already processed correctly, continue to next expt 
    if exist(TimingFile,'file') 
            load(TimingFile)
    else 
                return
    end 
    
    ThorSyncFile = fullfile(LocalPath, 'Episode001.h5');
    
    TempFile = 'C:\Users\kelsonss\Desktop\Temp\Episode001.h5';
    try
    copyfile(ThorSyncFile,TempFile)
    licks_raw=h5read(TempFile,'/AI/Licks');
    delete(TempFile)
    catch 
        return
    end 
    
    %licks = findpeaks(licks_raw,3);
    licks = findpeaks(diff(licks_raw),.4);
     
    licks = unique(floor(licks.loc/100) * 100 ) / 1000; 

    Out.RawLickFrames = licks;
    
    Out.TrialLickFrames = {};
   
    for trial_idx = 1:size(TimingInfo.FrameIdx,1)
        % get start and end of each trial
        t_start = TimingInfo.FrameIdx(trial_idx,1);
        t_end = TimingInfo.FrameIdx(trial_idx,2);
        % find licks tthat occured duing said trial
        lick_idx = (licks >= t_start)  & (licks < t_end);
        % shift lick times to be relative to tone onset
        t_licks =  licks(lick_idx) - t_start ;
        if isempty(t_licks)
            t_licks = nan;
        end 
        Out.TrialLickFrames{trial_idx} = t_licks;
    end 
        
        
  
    
 %% testing 
    figure
    subplot(1,4,1)
    all_expt_licks = cell2mat(cellfun( @(x) x(:)',Out.TrialLickFrames,'UniformOutput',0)) 
    histogram(all_expt_licks,[0:1:150])
 
  
   first_response = handles.FirstResponse(:,1);
   first_response_detect = cellfun(@(x) x(1),Out.TrialLickFrames)'/handles.pfs;
   
   [first_response, first_response_detect]
    
   min_diff = []
  for ii = 1:length(first_response)
    min_diff(ii) =  first_response(ii) -    
      
  
    %% saving 
    OutPath = fullfile(SavePath,'BehavioralResponses_Frames.mat');
    save(OutPath,'Out')
    fprintf('%d Minutes \n', floor(toc /60)) 
    
    

    

        

