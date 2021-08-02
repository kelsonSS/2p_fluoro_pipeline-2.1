function TimingInfo = getTimingInfo_H5(ThorSyncFile, PsignalTimingFile, fps,xml,num_frames_actual,debug)
%This function finds the frames that correspond to trials by looking
% comparing the ThorSync and Psignal Files)
% 

if ~exist('debug','var'); debug = 0;end 
if ~exist('fps','var');fps = 30 ; end
     
TimingInfo = struct()
%% Load Psignal timing file and Relevant Info
PsignalData = load(PsignalTimingFile);

TrialObjectClass = PsignalData.exptparams.TrialObjectClass;
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
first_frame = find(fc,1);
catch
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
             fprintf('\n %s is bad \n', ThorSyncFile)
             TimingInfo.Errors =  sprintf('\n %s is bad \n', ThorSyncFile);
             return
        end
        
    end
end



nchannels = xml.format{2}(3);



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
off = off(off > 3000); % ensure the first off trigger happens at least .1 sec after trial start



if length(on) ~=TotalTrials || length(on) ~= length(off)
    delta_peak = .1; % how much to change peak val by 
    pv_lower_flag= 0;  % if the last trial was too small
    pv_higher_flag = 0; % if last trial was too large
end 

while length(on) ~= TotalTrials || length(on) ~= length(off)
    m = length(on);
    n = length(off);
    g = length(gate);
    if m > 5000 || peak_val < 0.001 || peak_val > 2.5 || delta_peak < 5e-2
        warning('trial error manual inspection needed \n %s \n',ThorSyncFile)
        figure; plot(gate); hold on ; plot(fo);
        title(ThorSyncFile,'Interpreter','None')
        drawnow()
 %       if debug
           data_quality= questdlg('Is Data Usable if shortened?')
           if strmatch(data_quality,'Yes')  
              continue_plotting = 'Yes';
               while strmatch(continue_plotting, 'Yes'); 
               prompt = {'Enter expected peak size'};
               dlgtitle = 'Input';
               dims= [1 35];
               definput = {'.5'};
               prompt_peak_val = str2num(char(...
                             inputdlg(prompt,dlgtitle,dims,definput)));
               on =  findpeaks(diff(gate),prompt_peak_val);
               off = findpeaks(diff(gate * -1),prompt_peak_val);
               on = on.loc;
               off = off.loc;
               off = off(off > 3000); % ensure the first off trigger happens at least .1 sec after trial start
               m = length(on);
               figure; plot(diff(gate)); hold on ; %plot(diff(fo));
               x_lims = xlim;
               y_lims = ylim;
               plot([x_lims(1) x_lims(2)], [prompt_peak_val prompt_peak_val])
               scatter(on,repmat(y_lims(2),length(on),1)./2,'ko')
               scatter(off,repmat(y_lims(2),length(off),1)./2,'kx')
               fprintf('estimated trials:%d actual: %d peak val: %.2f  \n',...
            m,TotalTrials,prompt_peak_val)
               pause
               continue_plotting = questdlg('Try another value?');
               end 
               break

                         
               
           else 
               TimingInfo.Errors = 'Experiment was performed improperly'
               return
           end 
%        end 
          %  [SyncPath SyncName] = fileparts(ThorSyncFile);
          % SyncName = strcat(SyncName,'.h5');
          % LoadSyncEpisodeEdit(SyncPath, SyncName );
       %     pause 
      %  end 
       
    end
        
    
    
    % if we have switched directions decrease the delta size
    if pv_lower_flag && pv_higher_flag
        delta_peak = delta_peak/2;
        pv_lower_flag = 0;
        pv_higher_flag = 0;
    end 
    if m < TotalTrials || n < TotalTrials
        pv_lower_flag = 1;
        peak_val= peak_val - delta_peak;
        fprintf('too few estimated trials \n')   
        fprintf('estimated trials:%d actual: %d peak val: %.2f Delta: %.2f \n',...
            m,TotalTrials,peak_val,delta_peak)
    end 
    
      if m > TotalTrials || n> TotalTrials
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
    % ensure first off stimulus occurs at least a 10th of a second
    % after experiment start 
    off_idx = off_idx & (off' > 3000); 
    off = off(off_idx);
   
    
    
end  
    

frames = findpeaks(diff(fo),1);
frames = frames.loc;

trial_frames = cell(length(on),1);
for trial = 1:length(on)
    % find frames that happened after trial onset and before trial offset
    trial_frames_temp =  find( frames > on(trial) & frames<off(trial));
    % find the closest frame to the start of the trial 
   [~, closest_frame] = min(...
                        abs(...
                            frames-on(trial) ));
   
    
    % the closest_frame should either be the first frame or one before it 
    try
    assert( (closest_frame == trial_frames_temp(1)) ||...
             abs(closest_frame -trial_frames_temp(1))<=1)
    catch 
    end 
         
         
   % create list of frames in each trial
   if length(trial_frames_temp)>0
   if closest_frame ~= trial_frames_temp(1)
       trial_frames{trial} = cat(1,closest_frame,trial_frames_temp(1:end-1));
   else 
       trial_frames{trial} = trial_frames_temp;
   end 
   end
   % ensure consistent number of frames per trial 
   if length(trial_frames{trial})~= TimingInfo.tarFnums
       % too many frames
       if length(trial_frames{trial})> TimingInfo.tarFnums
           trial_frames{trial} = trial_frames{trial}(1:TimingInfo.tarFnums);
      % not enough frames 
      % this shouldn't happen unless an error
       else 
           warning('not enough frames found in trial %d',trial)
           % known cases could be
           %  1 frame missing on first trial (add one frame to beggining of movie and thorsync file) 
           %  user error closing early on last trial ( discard last trial)
           %  computer error midway through trial (discard expt)
       end 
   end 
   
end     

%% sanity checks

if num_frames_actual < TotalTrials * TimingInfo.tarFnums -1
    warning('trial has too few frames, skipping')
    return
end 


num_frames_predicted = findpeaks(diff(fo),1);
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


 frames_actual = cellfun(@length,trial_frames);
 bad_frame_flag = mode(frames_actual)  ~= TimingInfo.tarFnums ;
 
 bad_timing_flag = (floor(mean(off-on)/1000) ~= TimingInfo.tarFnums);
 
 % sanity check
 if bad_timing_flag || bad_frame_flag
     warning(sprintf('%s may be wrong file! manual inspection needed',ThorSyncFile))
     TimingInfo.Error = 'Suspected File Mismatch'
 end 


%% Create TimingInfo Structure

 

 
%Find frame timing for each trial; length(SeqEndVals)==#trials
FrameIdx(:,1) = cellfun(@min,trial_frames);
FrameIdx(:,2) = cellfun(@max,trial_frames) ;

TimingInfo.FrameIdx = FrameIdx;
TimingInfo.SeqEndVals = FrameIdx(:,2);
%TimingInfo.ITI = ITI;


% function [on off] = extract_h5_stimulus_mode(on,off)
% 
% 
% 
% if num_frames_actual == TimingInfo.tarFnums  * TotalTrials
%     on = 1 + TimingInfo.tarFnums .*(0:TotalTrials-1);
%     off = on + TimingInfo.tarFnums - 1; 
% return
% end
%     
% 
% end
% 
% 




        


