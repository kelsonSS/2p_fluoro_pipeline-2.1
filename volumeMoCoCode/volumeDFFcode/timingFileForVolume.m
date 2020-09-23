%% add timing file for volumetric datasets
clear

directory = 'C:\Users\Zac\Dropbox (Personal)\research\volume_data\031120_352ll_100um20st_FRA_diffxy\';

%% psignal info

% cd(directory)
% psigfile = dir('RND*.mat');
% handles =  WF_getPsignalInfo(fullfile(psigfile.folder,psigfile.name));

%% code from to look at which frame trigger crosses closest to another thorsync signal
% in this case it looks at the PsignalGate for trial start time

pathname = fullfile(directory,'sync001\');
filename = 'Episode001.h5';

%cd(pathname)
LoadSyncEpisodeEdit(pathname,filename)

if ~exist('ai2','var')
%     ai2 = Frame_In;
    ai2 = PsignalGate;
end
ai2 = ai2/max(ai2);
Frame_Out = Frame_Out/max(Frame_Out);
th = 0.95; % the threshold to detect pulses
ai2_crossing = find(ai2(1:end-1)<th & ai2(2:end)>=th);
Frame_Out_crossing = find(Frame_Out(1:end-1)<th & Frame_Out(2:end)>=th);
% for each ai2 crossing, find the closest Frame_Out
trialStarts = zeros(length(ai2_crossing),1);
for j = 1:length(ai2_crossing)
    dist = (Frame_Out_crossing - ai2_crossing(j)).^2;
    [~,trialStarts(j)] = min(dist);
end
trialStartFrame = trialStarts;
%         stimOnFrame = updateFrameTimes(pathname,filename);

cd ..
save(['trialStartFrames.mat'],'trialStartFrame')

