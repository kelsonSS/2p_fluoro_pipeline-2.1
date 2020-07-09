function CheckExperimentQuality(input)
%Performs an automated number of quality checks of the data and exports them 
% as summary figures to the analysis folder 


tic();
%Find index of experiment to be used for registration template
%regidx=find(strcmpi(input.expname,input.regexp));
%Path of images used for registration template
for expnum = 1:length(input.expname)

bb=strsplit(input.path,filesep);
newpath = char(strcat(input.savepath, filesep, bb{end} , filesep , input.expname{expnum}, filesep )) ;
disp(newpath)

if ~exist(fullfile(newpath,'Fluorescence.mat'),'file')
    continue
end 
load(fullfile(newpath,'Fluorescence.mat')) % Fluorescence
load(fullfile(newpath,'CellDefinitions.mat'))%  Cell Definitions
load(fullfile(newpath,'TimingInfo.mat'))% extracted Timing info

ThorSyncFile = fullfile(newpath,'Episode001.h5');
 %% load thorSync Data    
if exist(ThorSyncFile,'file')   
    % load Frameout
    try
    Sync.fo=h5read(ThorSyncFile,'/DI/Frame Out');

    catch    
  fprintf('%s is bad /n', ThorSyncFile)
    end 

% load counter if it was used
try
Sync.fc=h5read(ThorSyncFile,'/CI/Frame Counter'); % may be unreliable
catch
    Sync.fc = Sync.fo;
end 
  % load Behavioral Sync.gate ( high when trial is occuring)    
try
    Sync.gate=h5read(ThorSyncFile,'/AI/ai5'); % trial Sync.gate signal
catch
    try
        Sync.gate=h5read(ThorSyncFile,'/AI/PsignalSync.gate');
    catch
        try
        Sync.gate=h5read(ThorSyncFile,'/AI/ai4');
        catch
             fprintf('%s is bad /n', ThorSyncFile)
        end
        
    end
end
else 
    continue
end


% load XML
opts = get_options_from_xml(fullfile(newpath,'Experiment.xml'));

% load first minute of images 
fh = fopen(fullfile(newpath,'greenchannelregistered.raw'));

NumFrames = 8000;
ItemsPerImage = opts.dimX *opts.dimY ;
numFrames_total =  FindRawImgSize(fh,[opts.dimX opts.dimY]);
NumFrames = min(numFrames_total,NumFrames);
IMG =fread(fh,ItemsPerImage * NumFrames,'uint16');
fclose(fh);
IMG = reshape(IMG,opts.dimX,opts.dimY,[]);

% load Psignal info
handles = WF_getPsignalInfo(fullfile(newpath,input.psignalfiles{expnum}));

   on =  findpeaks(diff(Sync.gate),1);
      % ensure on peaks are acutally on

    on = on.loc;
    
 
    off = findpeaks(diff(Sync.gate * -1),1);
    %ensure off peaks are off  
    off = off.loc;

    
    
 frames_actual = findpeaks(diff(Sync.fo(on(1):off(1))),1);
 bad_frame_flag = length(frames_actual.loc) ~= TimingInfo.tarFnums ;
 
 bad_timing_flag = floor(mean(off-on)/1000) ~= TimingInfo.tarFnums;
 
    

 clear on 
 clear off
 clear off_idx
 clear on_idx

% plot Gate and FrameTiming for first 10 trials 

figure
plot(Sync.gate)
hold on 
plot(Sync.fo)


key_frames(1) = handles.PreStimSilence;
key_frames(2) = key_frames(1) +handles.PrimaryDuration;
key_frames(3) = key_frames(2) + handles.PostStimSilence;
key_names = {'Silence','Sound-On','Silence'};

if handles.BackgroundNoise(1) ~= -99 
     NoiseOn  = handles.BackgroundNoise(2) ;
     NoiseOFF = handles.BackgroundNoise(3)
% add more to addnoise to key_frames and add_names to Noise-on and noise
% off 

end



PrestimSilence = handles.pfs * handles.PreStimSilence;
SoundOn = PrestimSilence+  handles.pfs + handles.PrimaryDuration;

expt = 1;
Fpad = 30; 
figure('Renderer', 'painters', 'Position', [400 100 900 900])
while TimingInfo.FrameIdx(expt,2) < NumFrames
 
    % get the key_times for this experiment

    Fstart = TimingInfo.FrameIdx(expt,1);
    Fend  = TimingInfo.FrameIdx(expt,2);
     Fstart_actual = Fstart;
     Fend_actual = Fend;
     key_frames_expt = Fstart + key_frames .* handles.pfs;
     
    % play 1 second before and after each trial to ensure proper trial
    % aquisition
      
    
    Fstart = max( Fstart - Fpad, 1);
    Fend  = min( Fend + Fpad ,NumFrames);
    % find how each frame relates to thor sync timings
    frame_times = findpeaks(diff(Sync.fo),1);
    frame_times = frame_times.loc;
    
    
    
    for ii = Fstart:Fend
        trial_event = KeyTitle(ii,key_frames_expt,key_names);
        subplot(5,5,[1:20])
        imagesc(IMG(:,:,ii))
        colormap gray
        axis off 
        axis tight
        title(sprintf('Frame %d: %s ',ii-Fstart+1, trial_event),...
              'interpreter','none')
        
       subplot(5,1,5)
       plot(Sync.gate(frame_times(Fstart):frame_times(ii))) 
       pause(.03)
     
        
    end 
    cla



expt = expt+1; 
end
end

function title = KeyTitle(frame,key_times,key_names)

% key_times and names are key_value pairs  Key title comapres the current 
% frame to the key times and displays the first name with the highest value
% thats  less_than or equal to the key time 

 id = find(frame<= key_times,1,'first');
try
 title =  key_names{id} ;
catch
    title =  'trial over';
end 





