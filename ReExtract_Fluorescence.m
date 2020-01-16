function Output =  ReExtract_Fluorescence(data_dir)

if isstruct(data_dir)
    data_dir = data_dir.Experiment_files;
end 




for expt = 1:length(data_dir)
    % file loading & unpacking
    load(fullfile(data_dir{expt},'Fluorescence.mat'))
    load(fullfile(data_dir{expt},'TimingInfo.mat'))
    
    numframes = TimingInfo.tarFnums;
    reps = length(TimingInfo.SeqEndVals);
    neurons = size(Output.fluoAllRaw,2);
    on = TimingInfo.FrameIdx(:,1);
    off = TimingInfo.FrameIdx(:,2);
    
    
    % pre-allocation
    FCell = zeros(numframes,reps,neurons);
    FCellCorrected = zeros(numframes,reps,neurons);
    FNeuropil = zeros(numframes,reps,neurons);
    
    rec_frames =  off(1) - on(1) + 1; 
    
    if rec_frames ~= numframes
        frame_diff = rec_frames - numframes;
        off = off - frame_diff;
    end 
    
    % loading
    for rep = 1:reps
       FCell(:,rep,:) =      Output.fluoAllRaw(on(rep):off(rep),:);
       FCellCorrected(:,rep,:)= Output.NPfluoAll(on(rep):off(rep),:);   
       FNeuropil(:,rep,:)=  Output.fluoAllCorr(on(rep):off(rep),:);
    end
    
    
B_Vec = mean(Output.fluoAllCorr(1:30,:,:));
Vec_DFF = (Output.fluoAllCorr -B_Vec)./B_Vec * 100;   
bad = min(Vec_DFF) < -50 ; 
Vec_DFF(:,bad) = [];

figure
hold on 

title(data_dir(expt))
plot(Vec_DFF)

for  rep = 1:length(on)
plot([on(rep) off(rep)],[200 200])
end

figure 
plot(squeeze(nanmean(FCellCorrected,2)))
title( data_dir(expt))
end 


    
    
    %saving 
%     Output.FCell = Fcell;
%     Output.FCellCorrected = FCellCorrected;
%     Output.FNeuropil = FNeuropil;
%     
%     save('Fluorescence.mat','Output')
    
    
    
    
end 
    
    
